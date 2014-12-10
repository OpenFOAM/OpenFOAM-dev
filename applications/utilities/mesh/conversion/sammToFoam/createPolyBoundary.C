/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Create intermediate mesh files from SAMM files

\*---------------------------------------------------------------------------*/

#include "sammMesh.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void sammMesh::createPolyBoundary()
{
    label nBoundaryFacesFound = 0;

    polyBoundaryPatchStartIndices_.setSize(boundary_.size());

    label nCreatedFaces = nInternalFaces_;

    const labelListList& PointCells = pointCells();

    forAll(boundary_, patchI)
    {
        const faceList& curShapePatch = boundary_[patchI];

        polyBoundaryPatchStartIndices_[patchI] = nCreatedFaces;

        forAll(curShapePatch, faceI)
        {
            bool found = false;

            const face& curFace = curShapePatch[faceI];

            meshFaces_[nCreatedFaces] = curFace;

            // Must find which cell this face belongs to in order to
            // mark it in the cellPolys_
            const labelList& facePoints = curFace;

            forAll(facePoints, pointI)
            {
                const labelList& facePointCells =
                    PointCells[facePoints[pointI]];

                forAll(facePointCells, cellI)
                {
                    const faceList& curCellFaces =
                        cellFaces_[facePointCells[cellI]];

                    forAll(curCellFaces, cellFaceI)
                    {
                        if (curCellFaces[cellFaceI] == curFace)
                        {
                            // Found the cell face corresponding to this face
                            found = true;

                            // Debugging
                            if
                            (
                                cellPolys_[facePointCells[cellI]][cellFaceI]
                             != -1
                            )
                            {
                                FatalErrorIn
                                (
                                    "void sammMesh::createPolyBoundary()"
                                )   << "This looks like an already detected "
                                    << "internal face"
                                    << abort(FatalError);
                            }

                            cellPolys_[facePointCells[cellI]][cellFaceI] =
                                nCreatedFaces;

                            nBoundaryFacesFound++;
                        }
                        if (found) break;
                    }
                    if (found) break;
                }
                if (found) break;
            }

            nCreatedFaces++;
        }
    }

    // reset the size of the face list
    meshFaces_.setSize(nCreatedFaces);

    Info<< "Number of boundary faces: " << nBoundaryFacesFound << endl;
    Info<< "Total number of faces: " << nCreatedFaces << endl;
}


List<polyPatch* > sammMesh::polyBoundaryPatches(const polyMesh& pMesh)
{
    List<polyPatch* > p(boundary_.size());

    forAll(boundary_, patchI)
    {
        const faceList& curShapePatch = boundary_[patchI];

        p[patchI] = polyPatch::New
        (
            patchTypes_[patchI],
            patchNames_[patchI],
            curShapePatch.size(),
            polyBoundaryPatchStartIndices_[patchI],
            patchI,
            pMesh.boundaryMesh()
        ).ptr();

        p[patchI]->physicalType() = patchPhysicalTypes_[patchI];
    }

    return p;
}


// ************************************************************************* //
