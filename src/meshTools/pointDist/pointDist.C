/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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
#include "pointEdgeDist.H"
#include "pointMesh.H"
#include "PointEdgeWave.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointDist::pointDist
(
    const pointMesh& pMesh,
    const labelHashSet& patchIDs,
    const pointField& points,
    const scalar maxDist
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
    patchIndices_(patchIDs),
    maxDist_(maxDist)
{
    correct();
}


Foam::pointDist::pointDist
(
    const pointMesh& pMesh,
    const labelHashSet& patchIDs,
    const labelHashSet& zoneIDs,
    const pointField& points,
    const scalar maxDist
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
    patchIndices_(patchIDs),
    zoneIndices_(zoneIDs),
    maxDist_(maxDist)
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

    forAllConstIter(labelHashSet, patchIndices_, iter)
    {
        const label patchi = iter.key();
        nPatchPoints += pbm[patchi].meshPoints().size();
    }

    forAllConstIter(labelHashSet, zoneIndices_, iter)
    {
        const label zonei = iter.key();
        nPatchPoints += mesh()().pointZones()[zonei].size();
    }

    pointEdgeDist::data pointEdgeData(points_, maxDist_);

    // Set initial changed points to all the patch points(if patch present)
    List<pointEdgeDist> patchPointsInfo(nPatchPoints);
    labelList patchPoints(nPatchPoints);
    nPatchPoints = 0;

    // Add the patch points to the patchPointsInfo
    forAllConstIter(labelHashSet, patchIndices_, iter)
    {
        const label patchi = iter.key();
        const labelList& mp = pbm[patchi].meshPoints();

        forAll(mp, ppi)
        {
            const label meshPointi = mp[ppi];
            patchPoints[nPatchPoints] = meshPointi;
            patchPointsInfo[nPatchPoints] = pointEdgeDist
            (
                pointEdgeData.points[meshPointi],
                0
            );
            nPatchPoints++;
        }
    }

    // Add the zone points to the patchPointsInfo
    forAllConstIter(labelHashSet, zoneIndices_, iter)
    {
        const label zonei = iter.key();
        const labelList& zonePoints = mesh()().pointZones()[zonei];

        forAll(zonePoints, j)
        {
            const label meshPointi = zonePoints[j];
            patchPoints[nPatchPoints] = meshPointi;
            patchPointsInfo[nPatchPoints] = pointEdgeDist
            (
                pointEdgeData.points[meshPointi],
                0
            );
            nPatchPoints++;
        }
    }

    // Current info on points
    List<pointEdgeDist> allPointInfo(mesh()().nPoints());

    // Current info on edges
    List<pointEdgeDist> allEdgeInfo(mesh()().nEdges());

    PointEdgeWave
    <
        pointEdgeDist,
        pointEdgeDist::data
    > patchCalc
    (
        mesh()(),
        patchPoints,
        patchPointsInfo,

        allPointInfo,
        allEdgeInfo,
        mesh().globalData().nTotalPoints(), // max iterations
        pointEdgeData
    );

    pointScalarField& psf = *this;

    forAll(allPointInfo, pointi)
    {
        if (allPointInfo[pointi].valid(pointEdgeData))
        {
            psf[pointi] = sqrt(allPointInfo[pointi].distSqr());
        }
        else
        {
            psf[pointi] = maxDist_;
        }
    }

    forAllConstIter(labelHashSet, zoneIndices_, iter)
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
