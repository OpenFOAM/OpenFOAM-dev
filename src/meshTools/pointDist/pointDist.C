/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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
    const labelHashSet& startPatchIDs,
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
    startPatchIDs_(startPatchIDs),
    maxDist_(maxDist)
{
    correct();
}


Foam::pointDist::pointDist
(
    const pointMesh& pMesh,
    const labelHashSet& startPatchIDs,
    const labelHashSet& startZoneIDs,
    const labelHashSet& endPatchIDs,
    const labelHashSet& endZoneIDs,
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
    startPatchIDs_(startPatchIDs),
    startZoneIDs_(startZoneIDs),
    endPatchIDs_(endPatchIDs),
    endZoneIDs_(endZoneIDs),
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

    label nKnownPoints = 0;

    forAllConstIter(labelHashSet, startPatchIDs_, iter)
    {
        nKnownPoints += pbm[iter.key()].meshPoints().size();
    }

    forAllConstIter(labelHashSet, startZoneIDs_, iter)
    {
        nKnownPoints += mesh()().pointZones()[iter.key()].size();
    }

    forAllConstIter(labelHashSet, endPatchIDs_, iter)
    {
        nKnownPoints += pbm[iter.key()].meshPoints().size();
    }

    forAllConstIter(labelHashSet, endZoneIDs_, iter)
    {
        nKnownPoints += mesh()().pointZones()[iter.key()].size();
    }

    pointEdgeDist::data pointEdgeData(points_, maxDist_);

    // Set the start points of the wave
    List<pointEdgeDist> knownPointsInfo(nKnownPoints);
    labelList knownPoints(nKnownPoints);
    nKnownPoints = 0;

    // Add the patch start points to the knownPointsInfo
    forAllConstIter(labelHashSet, startPatchIDs_, iter)
    {
        const label patchi = iter.key();
        const labelList& mp = pbm[patchi].meshPoints();

        forAll(mp, ppi)
        {
            const label meshPointi = mp[ppi];
            knownPoints[nKnownPoints] = meshPointi;
            knownPointsInfo[nKnownPoints] = pointEdgeDist
            (
                pointEdgeData.points[meshPointi],
                0
            );
            nKnownPoints++;
        }
    }

    // Add the zone start points to the knownPointsInfo
    forAllConstIter(labelHashSet, startZoneIDs_, iter)
    {
        const label zonei = iter.key();
        const labelList& zonePoints = mesh()().pointZones()[zonei];

        forAll(zonePoints, j)
        {
            const label meshPointi = zonePoints[j];
            knownPoints[nKnownPoints] = meshPointi;
            knownPointsInfo[nKnownPoints] = pointEdgeDist
            (
                pointEdgeData.points[meshPointi],
                0
            );
            nKnownPoints++;
        }
    }

    // Add the patch end points to the knownPointsInfo
    forAllConstIter(labelHashSet, endPatchIDs_, iter)
    {
        const label patchi = iter.key();
        const labelList& mp = pbm[patchi].meshPoints();

        forAll(mp, ppi)
        {
            const label meshPointi = mp[ppi];
            knownPoints[nKnownPoints] = meshPointi;
            knownPointsInfo[nKnownPoints] = pointEdgeDist
            (
                pointEdgeData.points[meshPointi],
                -great
            );
            nKnownPoints++;
        }
    }

    // Add the zone end points to the knownPointsInfo
    forAllConstIter(labelHashSet, endZoneIDs_, iter)
    {
        const label zonei = iter.key();
        const labelList& zonePoints = mesh()().pointZones()[zonei];

        forAll(zonePoints, j)
        {
            const label meshPointi = zonePoints[j];
            knownPoints[nKnownPoints] = meshPointi;
            knownPointsInfo[nKnownPoints] = pointEdgeDist
            (
                pointEdgeData.points[meshPointi],
                -great
            );
            nKnownPoints++;
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
        knownPoints,
        knownPointsInfo,

        allPointInfo,
        allEdgeInfo,
        mesh().globalData().nTotalPoints(), // max iterations
        pointEdgeData
    );

    pointScalarField& psf = *this;

    forAll(allPointInfo, pointi)
    {
        if (allPointInfo[pointi].set(pointEdgeData))
        {
            psf[pointi] = sqrt(allPointInfo[pointi].distSqr());
        }
        else
        {
            psf[pointi] = maxDist_;
        }
    }

    forAllConstIter(labelHashSet, startZoneIDs_, iter)
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
