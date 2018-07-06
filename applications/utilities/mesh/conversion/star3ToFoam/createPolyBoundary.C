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
    Create intermediate mesh files from PROSTAR files

\*---------------------------------------------------------------------------*/

#include "starMesh.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::starMesh::createPolyBoundary()
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
                                if
                                (
                                    cellPolys_[facePointCells[celli]][cellFacei]
                                  > nInternalFaces_
                                )
                                {
                                    InfoInFunction
                                        << "Problem with face: " << curFace
                                        << "\nProbably multiple definitions "
                                        << "of a single boundary face. " << endl
                                        << "Other boundary face: "
                                        << curCellFaces[cellFacei]
                                        << endl;

                                    Info<< "PROSTAR Command: vset,news,vlis";
                                    forAll(curCellFaces[cellFacei], spI)
                                    {
                                        // check if the point is given by STAR
                                        // or created locally
                                        if
                                        (
                                            curCellFaces[cellFacei][spI] > -1
                                         && curCellFaces[cellFacei][spI]
                                                < starPointID_.size()
                                        )
                                        {
                                            Info<< ","
                                                << starPointID_
                                                 [curCellFaces[cellFacei][spI]];
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
                                    InfoInFunction
                                        << "Problem with face: " << curFace
                                        << "\nProbably trying to define a "
                                        << "boundary face on a previously "
                                        << "matched internal face. " << endl
                                        << "Internal face: "
                                        << curCellFaces[cellFacei]
                                        << endl;

                                    Info<< "PROSTAR Command: vset,news,vlis";
                                    forAll(curCellFaces[cellFacei], spI)
                                    {
                                        // check if the point is given by STAR
                                        // or created locally
                                        if
                                        (
                                            curCellFaces[cellFacei][spI] > -1
                                         && curCellFaces[cellFacei][spI]
                                                < starPointID_.size()
                                        )
                                        {
                                            Info<< ","
                                                << starPointID_
                                                 [curCellFaces[cellFacei][spI]];
                                        }
                                        else
                                        {
                                            Info<< ",???";
                                        }
                                    }
                                    Info<< " $ bset,add,vset,all" << endl;

                                }
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

    // check all cellPolys_ to see if there are any missing faces
    label nMissingFaceFound = 0;

    forAll(cellPolys_, celli)
    {
        const labelList& curFaces = cellPolys_[celli];

        forAll(curFaces, facei)
        {
            if (curFaces[facei] < 0)
            {
                const face& missingFace = cellFaces_[celli][facei];

                InfoInFunction
                    << "Missing face found in cell " << celli
                    << ".\nType: " << cellShapes_[celli].model().name()
                    << ". STAR cell number: " << starCellID_[celli]
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

    forAll(cellPolys_, celli)
    {
        const labelList& curFaces = cellPolys_[celli];

        forAll(curFaces, facei)
        {
            markupFaces[curFaces[facei]]++;
        }
    }

    for (label i = nInternalFaces_; i < markupFaces.size(); i++)
    {
        markupFaces[i]++;
    }

    label nProblemFacesFound = 0;

    forAll(markupFaces, facei)
    {
        if (markupFaces[facei] != 2)
        {
            const face& problemFace = meshFaces_[facei];

            InfoInFunction
                << "Problem with face " << facei << ": addressed "
                << markupFaces[facei] << " times (should be 2!). Face: "
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


Foam::List<Foam::polyPatch*>
Foam::starMesh::polyBoundaryPatches(const polyMesh& pMesh)
{
    List<polyPatch*> p(boundary_.size());

    forAll(boundary_, patchi)
    {
        p[patchi] = polyPatch::New
        (
            patchTypes_[patchi],
            patchNames_[patchi],
            boundary_[patchi].size(),
            polyBoundaryPatchStartIndices_[patchi],
            patchi,
            pMesh.boundaryMesh()
        ).ptr();
    }

    return p;
}


// ************************************************************************* //
