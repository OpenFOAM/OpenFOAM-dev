/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "multiSolidBody_pointMeshMover.H"
#include "syncTools.H"
#include "polyTopoChangeMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(multiSolidBody, 0);
    addToRunTimeSelectionTable
    (
        pointMeshMover,
        multiSolidBody,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointMeshMovers::multiSolidBody::updateZonePointIndices()
{
    forAll(zoneIndices_, zonei)
    {
        const cellZone& zoneCells = poly().cellZones()[zoneIndices_[zonei]];

        boolList pointInZone(poly().nPoints(), false);

        forAll(zoneCells, zoneCelli)
        {
            const cell& c = poly().cells()[zoneCells[zoneCelli]];
            forAll(c, cfi)
            {
                const face& f = poly().faces()[c[cfi]];
                forAll(f, fpi)
                {
                    pointInZone[f[fpi]] = true;
                }
            }
        }

        syncTools::syncPointList(poly(), pointInZone, orEqOp<bool>(), false);

        DynamicList<label> zonePoints(poly().nPoints());
        forAll(pointInZone, celli)
        {
            if (pointInZone[celli])
            {
                zonePoints.append(celli);
            }
        }

        zonePoints_[zonei].transfer(zonePoints);
    }

    // Check that points are not in multiple zones
    labelList pointZone(poly().nPoints(), -1);
    forAll(zoneIndices_, zonei)
    {
        forAll(zonePoints_[zonei], zonePointi)
        {
            const label pointi = zonePoints_[zonei][zonePointi];
            pointZone[pointi] =
                pointZone[pointi] == -1 || pointZone[pointi] == zonei
              ? zonei
              : -2;
        }
    }
    syncTools::syncPointList
    (
        poly(),
        pointZone,
        [](label& x, const label& y) { x = x == -1 || x == y ? y : -2; },
        label(-1)
    );
    const label errorPointi = findIndex(pointZone, -2);
    if (errorPointi != -1)
    {
        FatalErrorInFunction
            << "Point " << errorPointi << " at " << poly().points()[errorPointi]
            << " is in multiple moving zones. This is not allowed."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::multiSolidBody::multiSolidBody
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementPoints0(mesh, dict, typeName)
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
    updateZonePointIndices();

    transforms_.setSize(zonei);
    forAll(zoneIndices_, zonei)
    {
        transforms_[zonei] = SBMFs_[zonei].transformation();
    }

    forAll(zoneIndices_, zonei)
    {
        Info<< "Applying solid-body motion " << SBMFs_[zonei].type()
            << " to cellZone " << mesh.cellZones()[zoneIndices_[zonei]].name()
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::multiSolidBody::~multiSolidBody()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::pointMeshMovers::multiSolidBody::newPoints()
{
    tmp<pointField> ttransformedPts(new pointField(poly().points()));
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


void Foam::pointMeshMovers::multiSolidBody::topoChange
(
    const polyTopoChangeMap& map
)
{
    updateZonePointIndices();

    // pointMesh already updates registered pointFields

    const pointField& points = poly().points();

    pointField newPoints0(points);

    forAll(zoneIndices_, zonei)
    {
        forAll(zonePoints_[zonei], zonePointi)
        {
            const label pointi = zonePoints_[zonei][zonePointi];

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
                    transforms_[zonei].invTransformPoint(points[pointi]);
            }
        }
    }

    twoDCorrectPoints(newPoints0);

    // Move into base class storage and mark as to-be-written
    points0_.transfer(newPoints0);
    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = poly().time().name();
}


void Foam::pointMeshMovers::multiSolidBody::distribute
(
    const polyDistributionMap& map
)
{
    displacementPoints0::distribute(map);

    updateZonePointIndices();
}


void Foam::pointMeshMovers::multiSolidBody::mapMesh(const polyMeshMap& map)
{
    displacementPoints0::mapMesh(map);

    updateZonePointIndices();

    pointField& points0 = this->points0();

    forAll(zoneIndices_, zonei)
    {
        forAll(zonePoints_[zonei], zonePointi)
        {
            const label pointi = zonePoints_[zonei][zonePointi];

            points0[pointi] =
                    transforms_[zonei].invTransformPoint(points0[pointi]);
        }
    }

    twoDCorrectPoints(points0);

    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = poly().time().name();
}


// ************************************************************************* //
