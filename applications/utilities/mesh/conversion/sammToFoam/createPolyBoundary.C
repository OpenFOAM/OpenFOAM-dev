/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

void Foam::sammMesh::createPolyBoundary()
{
    label nBoundaryFacesFound = 0;

    polyBoundaryPatchStartIndices_.setSize(boundary_.size());

    label nCreatedFaces = nInternalFaces_;

    const labelListList& PointCells = pointCells();

    forAll(boundary_, patchi)
    {
        const faceList& curShapePatch = boundary_[patchi];

        polyBoundaryPatchStartIndices_[patchi] = nCreatedFaces;

        forAll(curShapePatch, facei)
        {
            bool found = false;

            const face& curFace = curShapePatch[facei];

            meshFaces_[nCreatedFaces] = curFace;

            // Must find which cell this face belongs to in order to
            // mark it in the cellPolys_
            const labelList& facePoints = curFace;

            forAll(facePoints, pointi)
            {
                const labelList& facePointCells =
                    PointCells[facePoints[pointi]];

                forAll(facePointCells, celli)
                {
                    const faceList& curCellFaces =
                        cellFaces_[facePointCells[celli]];

                    forAll(curCellFaces, cellFacei)
                    {
                        if (curCellFaces[cellFacei] == curFace)
                        {
                            // Found the cell face corresponding to this face
                            found = true;

                            // Debugging
                            if
                            (
                                cellPolys_[facePointCells[celli]][cellFacei]
                             != -1
                            )
                            {
                                FatalErrorInFunction
                                    << "This looks like an already detected "
                                    << "internal face"
                                    << abort(FatalError);
                            }

                            cellPolys_[facePointCells[celli]][cellFacei] =
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


Foam::List<Foam::polyPatch* > Foam::sammMesh::polyBoundaryPatches
(
    const polyMesh& pMesh
)
{
    List<polyPatch* > p(boundary_.size());

    forAll(boundary_, patchi)
    {
        const faceList& curShapePatch = boundary_[patchi];

        p[patchi] = polyPatch::New
        (
            patchTypes_[patchi],
            patchNames_[patchi],
            curShapePatch.size(),
            polyBoundaryPatchStartIndices_[patchi],
            patchi,
            pMesh.boundaryMesh()
        ).ptr();

        p[patchi]->physicalType() = patchPhysicalTypes_[patchi];
    }

    return p;
}


// ************************************************************************* //
