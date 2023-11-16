/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "pointDist.H"
#include "externalPointEdgePoint.H"
#include "pointMesh.H"
#include "PointEdgeWave.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointDist::pointDist
(
    const pointMesh& pMesh,
    const labelHashSet& patchIDs,
    const pointField& points
)
:
    pointScalarField
    (
        IOobject
        (
            "pointDistance",
            pMesh.db().time().name(),
            pMesh.db()
        ),
        pMesh,
        dimensionedScalar(dimLength, great)
    ),
    points_(points),
    patchIDs_(patchIDs),
    nUnset_(0)
{
    correct();
}


Foam::pointDist::pointDist
(
    const pointMesh& pMesh,
    const labelHashSet& patchIDs,
    const labelHashSet& zoneIDs,
    const pointField& points
)
:
    pointScalarField
    (
        IOobject
        (
            "pointDistance",
            pMesh.db().time().name(),
            pMesh.db()
        ),
        pMesh,
        dimensionedScalar(dimLength, great)
    ),
    points_(points),
    patchIDs_(patchIDs),
    zoneIDs_(zoneIDs),
    nUnset_(0)
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointDist::~pointDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointDist::correct()
{
    const pointBoundaryMesh& pbm = mesh().boundary();

    label nPatchPoints = 0;

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        const label patchi = iter.key();
        nPatchPoints += pbm[patchi].meshPoints().size();
    }

    forAllConstIter(labelHashSet, zoneIDs_, iter)
    {
        const label zonei = iter.key();
        nPatchPoints += mesh()().pointZones()[zonei].size();
    }

    externalPointEdgePoint::trackingData td(points_);

    // Set initial changed points to all the patch points(if patch present)
    List<externalPointEdgePoint> patchPointsInfo(nPatchPoints);
    labelList patchPoints(nPatchPoints);
    nPatchPoints = 0;

    // Add the patch points to the patchPointsInfo
    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        const label patchi = iter.key();
        const labelList& mp = pbm[patchi].meshPoints();

        forAll(mp, ppi)
        {
            const label meshPointi = mp[ppi];
            patchPoints[nPatchPoints] = meshPointi;
            patchPointsInfo[nPatchPoints] = externalPointEdgePoint
            (
                td.points_[meshPointi],
                0
            );
            nPatchPoints++;
        }
    }

    // Add the zone points to the patchPointsInfo
    forAllConstIter(labelHashSet, zoneIDs_, iter)
    {
        const label zonei = iter.key();
        const labelList& zonePoints = mesh()().pointZones()[zonei];

        forAll(zonePoints, j)
        {
            const label meshPointi = zonePoints[j];
            patchPoints[nPatchPoints] = meshPointi;
            patchPointsInfo[nPatchPoints] = externalPointEdgePoint
            (
                td.points_[meshPointi],
                0
            );
            nPatchPoints++;
        }
    }

    // Current info on points
    List<externalPointEdgePoint> allPointInfo(mesh()().nPoints());

    // Current info on edges
    List<externalPointEdgePoint> allEdgeInfo(mesh()().nEdges());

    PointEdgeWave
    <
        externalPointEdgePoint,
        externalPointEdgePoint::trackingData
    > patchCalc
    (
        mesh()(),
        patchPoints,
        patchPointsInfo,

        allPointInfo,
        allEdgeInfo,
        mesh().globalData().nTotalPoints(), // max iterations
        td
    );

    pointScalarField& psf = *this;

    forAll(allPointInfo, pointi)
    {
        if (allPointInfo[pointi].valid(td))
        {
            psf[pointi] = sqrt(allPointInfo[pointi].distSqr());
        }
        else
        {
            nUnset_++;
        }
    }

    forAllConstIter(labelHashSet, zoneIDs_, iter)
    {
        const label zonei = iter.key();
        const labelList& zonePoints = mesh()().pointZones()[zonei];

        forAll(zonePoints, j)
        {
            const label meshPointi = zonePoints[j];
            psf[meshPointi] = 0;
        }
    }
}


// ************************************************************************* //
