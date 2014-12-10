/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "cellDistFuncs.H"
#include "polyMesh.H"
#include "wallPolyPatch.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellDistFuncs, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Find val in first nElems elements of elems.
Foam::label Foam::cellDistFuncs::findIndex
(
    const label nElems,
    const labelList& elems,
    const label val
)
{
    for (label i = 0; i < nElems; i++)
    {
        if (elems[i] == val)
        {
            return i;
        }
    }
    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellDistFuncs::cellDistFuncs(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelHashSet Foam::cellDistFuncs::getPatchIDs
(
    const wordReList& patchNames
) const
{
    return mesh().boundaryMesh().patchSet(patchNames, false);
}


// Return smallest true distance from p to any of wallFaces.
// Note that even if normal hits face we still check other faces.
// Note that wallFaces is untruncated and we explicitly pass in size.
Foam::scalar Foam::cellDistFuncs::smallestDist
(
    const point& p,
    const polyPatch& patch,
    const label nWallFaces,
    const labelList& wallFaces,
    label& minFaceI
) const
{
    const pointField& points = patch.points();

    scalar minDist = GREAT;
    minFaceI = -1;

    for (label wallFaceI = 0; wallFaceI < nWallFaces; wallFaceI++)
    {
        label patchFaceI = wallFaces[wallFaceI];

        pointHit curHit = patch[patchFaceI].nearestPoint(p, points);

        if (curHit.distance() < minDist)
        {
            minDist = curHit.distance();
            minFaceI = patch.start() + patchFaceI;
        }
    }

    return minDist;
}


// Get point neighbours of faceI (including faceI). Returns number of faces.
// Note: does not allocate storage but does use linear search to determine
// uniqueness. For polygonal faces this might be quite inefficient.
Foam::label Foam::cellDistFuncs::getPointNeighbours
(
    const primitivePatch& patch,
    const label patchFaceI,
    labelList& neighbours
) const
{
    label nNeighbours = 0;

    // Add myself
    neighbours[nNeighbours++] = patchFaceI;

    // Add all face neighbours
    const labelList& faceNeighbours = patch.faceFaces()[patchFaceI];

    forAll(faceNeighbours, faceNeighbourI)
    {
        neighbours[nNeighbours++] = faceNeighbours[faceNeighbourI];
    }

    // Remember part of neighbours that contains edge-connected faces.
    label nEdgeNbs = nNeighbours;


    // Add all point-only neighbours by linear searching in edge neighbours.
    // Assumes that point-only neighbours are not using multiple points on
    // face.

    const face& f = patch.localFaces()[patchFaceI];

    forAll(f, fp)
    {
        label pointI = f[fp];

        const labelList& pointNbs = patch.pointFaces()[pointI];

        forAll(pointNbs, nbI)
        {
            label faceI = pointNbs[nbI];

            // Check for faceI in edge-neighbours part of neighbours
            if (findIndex(nEdgeNbs, neighbours, faceI) == -1)
            {
                neighbours[nNeighbours++] = faceI;
            }
        }
    }


    if (debug)
    {
        // Check for duplicates

        // Use hashSet to determine nbs.
        labelHashSet nbs(4*f.size());

        forAll(f, fp)
        {
            const labelList& pointNbs = patch.pointFaces()[f[fp]];

            forAll(pointNbs, i)
            {
                nbs.insert(pointNbs[i]);
            }
        }

        // Subtract ours.
        for (label i = 0; i < nNeighbours; i++)
        {
            label nb = neighbours[i];

            if (!nbs.found(nb))
            {
                SeriousErrorIn("Foam::cellDistFuncs::getPointNeighbours")
                    << "getPointNeighbours : patchFaceI:" << patchFaceI
                    << " verts:" << f << endl;

                forAll(f, fp)
                {
                    SeriousErrorIn("Foam::cellDistFuncs::getPointNeighbours")
                        << "point:" << f[fp] << " pointFaces:"
                        << patch.pointFaces()[f[fp]] << endl;
                }

                for (label i = 0; i < nNeighbours; i++)
                {
                    SeriousErrorIn("Foam::cellDistFuncs::getPointNeighbours")
                        << "fast nbr:" << neighbours[i]
                        << endl;
                }

                FatalErrorIn("getPointNeighbours")
                    << "Problem: fast pointNeighbours routine included " << nb
                    << " which is not in proper neigbour list " << nbs.toc()
                    << abort(FatalError);
            }
            nbs.erase(nb);
        }

        if (nbs.size())
        {
            FatalErrorIn("getPointNeighbours")
                << "Problem: fast pointNeighbours routine did not find "
                << nbs.toc() << abort(FatalError);
        }
    }

    return nNeighbours;
}


// size of largest patch (out of supplied subset of patches)
Foam::label Foam::cellDistFuncs::maxPatchSize
(
    const labelHashSet& patchIDs
) const
{
    label maxSize = 0;

    forAll(mesh().boundaryMesh(), patchI)
    {
        if (patchIDs.found(patchI))
        {
            const polyPatch& patch = mesh().boundaryMesh()[patchI];

            maxSize = Foam::max(maxSize, patch.size());
        }
    }
    return maxSize;
}


// sum of patch sizes (out of supplied subset of patches)
Foam::label Foam::cellDistFuncs::sumPatchSize
(
    const labelHashSet& patchIDs
)
const
{
    label sum = 0;

    forAll(mesh().boundaryMesh(), patchI)
    {
        if (patchIDs.found(patchI))
        {
            const polyPatch& patch = mesh().boundaryMesh()[patchI];

            sum += patch.size();
        }
    }
    return sum;
}


// Gets nearest wall for cells next to wall
void Foam::cellDistFuncs::correctBoundaryFaceCells
(
    const labelHashSet& patchIDs,
    scalarField& wallDistCorrected,
    Map<label>& nearestFace
) const
{
    // Size neighbours array for maximum possible (= size of largest patch)
    label maxPointNeighbours = maxPatchSize(patchIDs);

    labelList neighbours(maxPointNeighbours);


    // Correct all cells with face on wall
    const vectorField& cellCentres = mesh().cellCentres();
    const labelList& faceOwner = mesh().faceOwner();

    forAll(mesh().boundaryMesh(), patchI)
    {
        if (patchIDs.found(patchI))
        {
            const polyPatch& patch = mesh().boundaryMesh()[patchI];

            // Check cells with face on wall
            forAll(patch, patchFaceI)
            {
                label nNeighbours = getPointNeighbours
                (
                    patch,
                    patchFaceI,
                    neighbours
                );

                label cellI = faceOwner[patch.start() + patchFaceI];

                label minFaceI = -1;

                wallDistCorrected[cellI] = smallestDist
                (
                    cellCentres[cellI],
                    patch,
                    nNeighbours,
                    neighbours,
                    minFaceI
                );

                // Store wallCell and its nearest neighbour
                nearestFace.insert(cellI, minFaceI);
            }
        }
    }

}


// Correct all cells connected to wall (via point) and not in nearestFace
void Foam::cellDistFuncs::correctBoundaryPointCells
(
    const labelHashSet& patchIDs,
    scalarField& wallDistCorrected,
    Map<label>& nearestFace
) const
{
    // Correct all (non-visited) cells with point on wall

    const vectorField& cellCentres = mesh().cellCentres();

    forAll(mesh().boundaryMesh(), patchI)
    {
        if (patchIDs.found(patchI))
        {
            const polyPatch& patch = mesh().boundaryMesh()[patchI];

            const labelList& meshPoints = patch.meshPoints();
            const labelListList& pointFaces = patch.pointFaces();

            forAll(meshPoints, meshPointI)
            {
                label vertI = meshPoints[meshPointI];

                const labelList& neighbours = mesh().pointCells(vertI);

                forAll(neighbours, neighbourI)
                {
                    label cellI = neighbours[neighbourI];

                    if (!nearestFace.found(cellI))
                    {
                        const labelList& wallFaces = pointFaces[meshPointI];

                        label minFaceI = -1;

                        wallDistCorrected[cellI] = smallestDist
                        (
                            cellCentres[cellI],
                            patch,
                            wallFaces.size(),
                            wallFaces,
                            minFaceI
                        );

                        // Store wallCell and its nearest neighbour
                        nearestFace.insert(cellI, minFaceI);
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
