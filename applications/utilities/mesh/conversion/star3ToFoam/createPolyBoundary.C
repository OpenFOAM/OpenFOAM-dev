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
    Create intermediate mesh files from PROSTAR files

\*---------------------------------------------------------------------------*/

#include "starMesh.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void starMesh::createPolyBoundary()
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
                                if
                                (
                                    cellPolys_[facePointCells[cellI]][cellFaceI]
                                  > nInternalFaces_
                                )
                                {
                                    Info
                                        << "void starMesh::createPolyBoundary()"
                                        << ": Problem with face: " << curFace
                                        << "\nProbably multiple definitions "
                                        << "of a single boundary face. " << endl
                                        << "Other boundary face: "
                                        << curCellFaces[cellFaceI]
                                        << endl;

                                    Info<< "PROSTAR Command: vset,news,vlis";
                                    forAll(curCellFaces[cellFaceI], spI)
                                    {
                                        // check if the point is given by STAR
                                        // or created locally
                                        if
                                        (
                                            curCellFaces[cellFaceI][spI] > -1
                                         && curCellFaces[cellFaceI][spI]
                                                < starPointID_.size()
                                        )
                                        {
                                            Info
                                                << ","
                                                << starPointID_
                                                 [curCellFaces[cellFaceI][spI]];
                                        }
                                        else
                                        {
                                            Info<< ",???";
                                        }
                                    }
                                    Info<< " $ bset,add,vset,all" << endl;
                                }
                                else
                                {
                                    Info
                                        << "void starMesh::createPolyBoundary()"
                                        << ": Problem with face: " << curFace
                                        << "\nProbably trying to define a "
                                        << "boundary face on a previously "
                                        << "matched internal face. " << endl
                                        << "Internal face: "
                                        << curCellFaces[cellFaceI]
                                        << endl;

                                    Info<< "PROSTAR Command: vset,news,vlis";
                                    forAll(curCellFaces[cellFaceI], spI)
                                    {
                                        // check if the point is given by STAR
                                        // or created locally
                                        if
                                        (
                                            curCellFaces[cellFaceI][spI] > -1
                                         && curCellFaces[cellFaceI][spI]
                                                < starPointID_.size()
                                        )
                                        {
                                            Info
                                                << ","
                                                << starPointID_
                                                 [curCellFaces[cellFaceI][spI]];
                                        }
                                        else
                                        {
                                            Info<< ",???";
                                        }
                                    }
                                    Info<< " $ bset,add,vset,all" << endl;

                                }
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

    // check all cellPolys_ to see if there are any missing faces
    label nMissingFaceFound = 0;

    forAll(cellPolys_, cellI)
    {
        const labelList& curFaces = cellPolys_[cellI];

        forAll(curFaces, faceI)
        {
            if (curFaces[faceI] < 0)
            {
                const face& missingFace = cellFaces_[cellI][faceI];

                Info<< "starMesh::createPolyBoundary() : "
                    << "missing face found in cell " << cellI
                    << ".\nType: " << cellShapes_[cellI].model().name()
                    << ". STAR cell number: " << starCellID_[cellI]
                    << ". Face: " << missingFace << endl;

                nMissingFaceFound++;

                Info<< "PROSTAR Command: vset,news,vlis";
                forAll(missingFace, spI)
                {
                    // check if the point is given by STAR or created locally
                    if
                    (
                        missingFace[spI] > -1
                     && missingFace[spI] < starPointID_.size()
                    )
                    {
                        Info<< "," << starPointID_[missingFace[spI]];
                    }
                    else
                    {
                        Info<< ",???";
                    }
                }
                Info<< " $ bset,add,vset,all" << endl;
            }
        }
    }

    if (nMissingFaceFound > 0)
    {
        Info<< "Number of unmatched faces: " << nMissingFaceFound << endl;
    }

    // reset the size of the face list
    meshFaces_.setSize(nCreatedFaces);

    // check the mesh for face mismatch
    // (faces addressed once or more than twice)
    labelList markupFaces(meshFaces_.size(), 0);

    forAll(cellPolys_, cellI)
    {
        const labelList& curFaces = cellPolys_[cellI];

        forAll(curFaces, faceI)
        {
            markupFaces[curFaces[faceI]]++;
        }
    }

    for (label i = nInternalFaces_; i < markupFaces.size(); i++)
    {
        markupFaces[i]++;
    }

    label nProblemFacesFound = 0;

    forAll(markupFaces, faceI)
    {
        if (markupFaces[faceI] != 2)
        {
            const face& problemFace = meshFaces_[faceI];

            Info<< "starMesh::createPolyBoundary() : "
                << "problem with face " << faceI << ": addressed "
                << markupFaces[faceI] << " times (should be 2!). Face: "
                << problemFace << endl;

            nProblemFacesFound++;

            Info<< "PROSTAR Command: vset,news,vlis";
            forAll(problemFace, spI)
            {
                // check if the point is given by STAR or created locally
                if
                (
                    problemFace[spI] > -1
                 && problemFace[spI] < starPointID_.size()
                )
                {
                    Info<< "," << starPointID_[problemFace[spI]];
                }
                else
                {
                    Info<< ",???";
                }
            }
            Info<< " $ bset,add,vset,all" << endl;
        }
    }

    if (nProblemFacesFound > 0)
    {
        Info<< "Number of incorrectly matched faces: "
            << nProblemFacesFound << endl;
    }

    Info<< "Number of boundary faces: " << nBoundaryFacesFound << endl;
    Info<< "Total number of faces: " << nCreatedFaces << endl;
}


List<polyPatch*> starMesh::polyBoundaryPatches(const polyMesh& pMesh)
{
    List<polyPatch*> p(boundary_.size());

    forAll(boundary_, patchI)
    {
        p[patchI] = polyPatch::New
        (
            patchTypes_[patchI],
            patchNames_[patchI],
            boundary_[patchI].size(),
            polyBoundaryPatchStartIndices_[patchI],
            patchI,
            pMesh.boundaryMesh()
        ).ptr();
    }

    return p;
}


// ************************************************************************* //
