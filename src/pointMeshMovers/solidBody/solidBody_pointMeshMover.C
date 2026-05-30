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

#include "solidBody_pointMeshMover.H"
#include "syncTools.H"
#include "polyTopoChangeMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(solidBody, 0);
    addToRunTimeSelectionTable
    (
        pointMeshMover,
        solidBody,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointMeshMovers::solidBody::updateZonePointIndices()
{
    zone_.regenerate();

    if (zone_.all())
    {
        zonePoints_.clear();
        return;
    }

    boolList pointInZone(poly().nPoints(), false);

    forAll(zone_.zone(), zoneCelli)
    {
        const cell& c = poly().cells()[zone_.zone()[zoneCelli]];
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

Foam::pointMeshMovers::solidBody::solidBody
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementPoints0(mesh, dict, typeName),
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

Foam::pointMeshMovers::solidBody::~solidBody()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::pointMeshMovers::solidBody::newPoints()
{
    transform_ = SBMFPtr_().transformation();

    if (zone_.all())
    {
        return transformPoints(transform_, points0_);
    }
    else
    {
        tmp<pointField> ttransformedPts(new pointField(poly().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, zonePoints_) = transformPoints
        (
            transform_,
            pointField(points0_, zonePoints_)
        );

        return ttransformedPts;
    }
}


void Foam::pointMeshMovers::solidBody::topoChange(const polyTopoChangeMap& map)
{
    zone_.topoChange(map);
    updateZonePointIndices();

    // pointMesh already updates pointFields

    const pointField& points = poly().points();

    pointField newPoints0(points);

    const label nZonePoints =
        zone_.all() ? poly().nPoints() : zonePoints_.size();

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
    points0_.instance() = poly().time().name();
}


void Foam::pointMeshMovers::solidBody::distribute
(
    const polyDistributionMap& map
)
{
    displacementPoints0::distribute(map);

    zone_.distribute(map);
    updateZonePointIndices();
}


void Foam::pointMeshMovers::solidBody::mapMesh(const polyMeshMap& map)
{
    displacementPoints0::mapMesh(map);

    zone_.mapMesh(map);
    updateZonePointIndices();

    pointField& points0 = this->points0();

    const label nZonePoints =
        zone_.all() ? poly().nPoints() : zonePoints_.size();

    for (label zonePointi = 0; zonePointi < nZonePoints; ++ zonePointi)
    {
        const label pointi =
            zone_.all() ? zonePointi : zonePoints_[zonePointi];

        points0[pointi] = transform_.invTransformPoint(points0[pointi]);
    }

    twoDCorrectPoints(points0);

    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = poly().time().name();
}


// ************************************************************************* //
