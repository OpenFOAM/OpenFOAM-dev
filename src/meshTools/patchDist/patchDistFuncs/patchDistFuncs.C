/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "patchDistFuncs.H"
#include "polyMesh.H"
#include "wallPolyPatch.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
Foam::patchDistFuncs::NoMap<Foam::labelPair>
Foam::patchDistFuncs::NoMap<Foam::labelPair>::null =
Foam::patchDistFuncs::NoMap<Foam::labelPair>();


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::patchDistFuncs::smallestDist
(
    const point& p,
    const polyPatch& patch,
    const labelUList& wallFaces,
    label& minPatchFacei
)
{
    const pointField& points = patch.points();

    scalar minDist = great;
    minPatchFacei = -1;

    forAll(wallFaces, wallFacei)
    {
        const label patchFacei = wallFaces[wallFacei];

        pointHit curHit = patch[patchFacei].nearestPoint(p, points);

        if (curHit.distance() < minDist)
        {
            minDist = curHit.distance();
            minPatchFacei = patchFacei;
        }
    }

    return minDist;
}


void Foam::patchDistFuncs::getPointNeighbours
(
    const primitivePatch& patch,
    const label patchFacei,
    DynamicList<label>& neighbours
)
{
    neighbours.clear();

    // Add myself
    neighbours.append(patchFacei);

    // Add all face neighbours
    const labelList& faceNeighbours = patch.faceFaces()[patchFacei];
    forAll(faceNeighbours, faceNeighbourI)
    {
        neighbours.append(faceNeighbours[faceNeighbourI]);
    }

    // Remember part of neighbours that contains edge-connected faces.
    const label nEdgeNbs = neighbours.size();

    // Add all point-only neighbours by linear searching in edge neighbours.
    // Assumes that point-only neighbours are not using multiple points on
    // face.
    const face& f = patch.localFaces()[patchFacei];
    forAll(f, fp)
    {
        const label pointi = f[fp];

        const labelList& pointNbs = patch.pointFaces()[pointi];

        forAll(pointNbs, nbI)
        {
            const label facei = pointNbs[nbI];

            // Check for facei in edge-neighbours part of neighbours
            if (findIndex(SubList<label>(neighbours, nEdgeNbs), facei) == -1)
            {
                neighbours.append(facei);
            }
        }
    }
}


Foam::label Foam::patchDistFuncs::maxPatchSize
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs
)
{
    label maxSize = 0;

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        maxSize = Foam::max(maxSize, mesh.boundaryMesh()[iter.key()].size());
    }

    return maxSize;
}


// ************************************************************************* //
