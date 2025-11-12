/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "multiSolidBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "cellZoneList.H"
#include "boolList.H"
#include "syncTools.H"
#include "polyTopoChangeMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiSolidBodyMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        multiSolidBodyMotionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiSolidBodyMotionSolver::updateZonePointIndices()
{
    forAll(zoneIndices_, zonei)
    {
        const cellZone& zoneCells = mesh().cellZones()[zoneIndices_[zonei]];

        boolList pointInZone(mesh().nPoints(), false);

        forAll(zoneCells, zoneCelli)
        {
            const cell& c = mesh().cells()[zoneCells[zoneCelli]];
            forAll(c, cfi)
            {
                const face& f = mesh().faces()[c[cfi]];
                forAll(f, fpi)
                {
                    pointInZone[f[fpi]] = true;
                }
            }
        }

        syncTools::syncPointList(mesh(), pointInZone, orEqOp<bool>(), false);

        DynamicList<label> zonePoints(mesh().nPoints());
        forAll(pointInZone, celli)
        {
            if (pointInZone[celli])
            {
                zonePoints.append(celli);
            }
        }

        zonePoints_[zonei].transfer(zonePoints);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionSolver::multiSolidBodyMotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    points0MotionSolver(name, mesh, dict, typeName)
{
    zoneIndices_.setSize(dict.size());
    SBMFs_.setSize(dict.size());
    label zonei = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (!iter().isDict()) continue;

        zoneIndices_[zonei] = mesh.cellZones().findIndex(iter().keyword());

        if (zoneIndices_[zonei] == -1)
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find cellZone named " << iter().keyword()
                << ". Valid zones are " << mesh.cellZones().toc()
                << exit(FatalIOError);
        }

        const dictionary& subDict = iter().dict();

        SBMFs_.set
        (
            zonei,
            solidBodyMotionFunction::New(subDict, mesh.time())
        );

        zonei ++;
    }
    zoneIndices_.setSize(zonei);
    SBMFs_.setSize(zonei);
    zonePoints_.setSize(zonei);
    transforms_.setSize(zonei);

    updateZonePointIndices();

    forAll(zoneIndices_, zonei)
    {
        Info<< "Applying solid body motion " << SBMFs_[zonei].type()
            << " to " << zonePoints_[zonei].size()
            << " points of cellZone "
            << mesh.cellZones()[zoneIndices_[zonei]].name() << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionSolver::~multiSolidBodyMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::multiSolidBodyMotionSolver::curPoints() const
{
    tmp<pointField> ttransformedPts(new pointField(mesh().points()));
    pointField& transformedPts = ttransformedPts.ref();

    forAll(zoneIndices_, zonei)
    {
        transforms_[zonei] = SBMFs_[zonei].transformation();

        UIndirectList<point>(transformedPts, zonePoints_[zonei]) =
            transformPoints
            (
                transforms_[zonei],
                pointField(points0_, zonePoints_[zonei])
            );
    }

    return ttransformedPts;
}


void Foam::multiSolidBodyMotionSolver::topoChange(const polyTopoChangeMap& map)
{
    updateZonePointIndices();

    const pointField& points = mesh().points();

    pointField newPoints0(map.pointMap().size());

    forAll(zoneIndices_, zonei)
    {
        boolList pointInZone(mesh().nPoints(), false);
        UIndirectList<bool>(pointInZone, zonePoints_[zonei]) = true;

        forAll(newPoints0, pointi)
        {
            const label oldPointi = map.pointMap()[pointi];

            if (oldPointi < 0)
            {
                FatalErrorInFunction
                    << "Cannot determine co-ordinates of introduced vertices."
                    << " New vertex " << pointi << " at co-ordinate "
                    << points[pointi] << exit(FatalError);
            }

            if (map.reversePointMap()[oldPointi] == pointi)
            {
                newPoints0[pointi] = points0_[oldPointi];
            }
            else
            {
                newPoints0[pointi] =
                    pointInZone[pointi]
                  ? transforms_[zonei].invTransformPoint(points[pointi])
                  : points[pointi];
            }
        }
    }

    twoDCorrectPoints(newPoints0);

    // Move into base class storage and mark as to-be-written
    points0_.transfer(newPoints0);
    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = mesh().time().name();
}


void Foam::multiSolidBodyMotionSolver::distribute
(
    const polyDistributionMap& map
)
{
    points0MotionSolver::distribute(map);

    updateZonePointIndices();
}


void Foam::multiSolidBodyMotionSolver::mapMesh(const polyMeshMap& map)
{
    points0MotionSolver::mapMesh(map);

    updateZonePointIndices();
}


// ************************************************************************* //
