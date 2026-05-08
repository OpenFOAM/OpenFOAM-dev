/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "solidBodyMeshMotion.H"
#include "generatedCellZone.H"
#include "transformField.H"
#include "syncTools.H"
#include "polyTopoChangeMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyMeshMotion, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        solidBodyMeshMotion,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidBodyMeshMotion::updateZonePointIndices()
{
    if (zone_.all())
    {
        zonePoints_.clear();
        return;
    }

    boolList pointInZone(mesh().nPoints(), false);

    forAll(zone_.zone(), zoneCelli)
    {
        const cell& c = mesh().cells()[zone_.zone()[zoneCelli]];
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

    zonePoints_.resize(count(pointInZone, true));

    label zonePointi = 0;
    forAll(pointInZone, pointi)
    {
        if (pointInZone[pointi])
        {
            zonePoints_[zonePointi ++] = pointi;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMeshMotion::solidBodyMeshMotion
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    points0MotionSolver(name, mesh, dict, typeName),
    SBMFPtr_(solidBodyMotionFunction::New(dict, mesh.time())),
    zone_(mesh, dict),
    zonePoints_(),
    transform_(SBMFPtr_().transformation())
{
    if (zone_.all())
    {
        Info<< "Applying solid body motion to entire mesh" << endl;
    }

    updateZonePointIndices();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMeshMotion::~solidBodyMeshMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::solidBodyMeshMotion::newPoints()
{
    transform_ = SBMFPtr_().transformation();

    if (zone_.all())
    {
        return transformPoints(transform_, points0_);
    }
    else
    {
        tmp<pointField> ttransformedPts(new pointField(mesh().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, zonePoints_) = transformPoints
        (
            transform_,
            pointField(points0_, zonePoints_)
        );

        return ttransformedPts;
    }
}


void Foam::solidBodyMeshMotion::topoChange(const polyTopoChangeMap& map)
{
    zone_.topoChange(map);
    updateZonePointIndices();

    // pointMesh already updates pointFields

    const pointField& points = mesh().points();

    pointField newPoints0(map.pointMap().size());

    const label nZonePoints =
        zone_.all() ? mesh().nPoints() : zonePoints_.size();

    for (label zonePointi = 0; zonePointi < nZonePoints; ++ zonePointi)
    {
        const label pointi =
            zone_.all() ? zonePointi : zonePoints_[zonePointi];

        const label oldPointi = map.pointMap()[pointi];

        if (oldPointi < 0)
        {
            FatalErrorInFunction
                << "Cannot determine co-ordinates of introduced points."
                << " New point " << pointi << " at " << points[pointi]
                << exit(FatalError);
        }

        if (map.reversePointMap()[oldPointi] == pointi)
        {
            newPoints0[pointi] = points0_[oldPointi];
        }
        else
        {
            newPoints0[pointi] =
                transform_.invTransformPoint(points[pointi]);
        }
    }

    twoDCorrectPoints(newPoints0);

    // Move into base class storage and mark as to-be-written
    points0_.transfer(newPoints0);
    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = mesh().time().name();
}


void Foam::solidBodyMeshMotion::distribute(const polyDistributionMap& map)
{
    points0MotionSolver::distribute(map);

    zone_.distribute(map);
    updateZonePointIndices();
}


void Foam::solidBodyMeshMotion::mapMesh(const polyMeshMap& map)
{
    points0MotionSolver::mapMesh(map);

    zone_.mapMesh(map);
    updateZonePointIndices();
}


// ************************************************************************* //
