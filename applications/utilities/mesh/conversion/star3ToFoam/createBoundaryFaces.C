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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialist version of face comparison to deal with
// PROSTAR boundary format idiosyncracies
bool starMesh::starEqualFace
(
    const face& boundaryFace,
    const face& cellFace
) const
{
    // A PROSTAR boundary face is defined by 4 vertices irrespective
    // of its topology.
    // In order to deal with all possibilities, cell face is
    // considered equal if three of the vertices are the same.
    bool cellFaceHappy = false;

    label nEqual = 0;

    forAll(cellFace, cellFaceLabelI)
    {
        const label curCellFaceLabel = cellFace[cellFaceLabelI];

        forAll(boundaryFace, bouFaceLabelI)
        {
            if (boundaryFace[bouFaceLabelI] == curCellFaceLabel)
            {
                nEqual++;

                break;
            }
        }
    }

    if (nEqual >= 3)
    {
        cellFaceHappy = true;
    }

    // Boundary face is happy if all of its vertices are recognised
    bool boundaryFaceHappy = true;

    forAll(boundaryFace, bouFaceLabelI)
    {
        const label curBouFaceLabel = boundaryFace[bouFaceLabelI];

        bool found = false;

        forAll(cellFace, cellFaceLabelI)
        {
            if (curBouFaceLabel == cellFace[cellFaceLabelI])
            {
                found = true;
                break;
            }
        }

        boundaryFaceHappy = boundaryFaceHappy && found;
    }

    return (cellFaceHappy && boundaryFaceHappy);
}


void starMesh::markBoundaryFaces()
{
    // set size of mark lists for the boundary
    boundaryCellIDs_.setSize(boundary_.size());
    boundaryCellFaceIDs_.setSize(boundary_.size());

    forAll(boundary_, patchI)
    {
        const faceList& patchFaces = boundary_[patchI];

        // set size of patch lists
        labelList& curBoundaryCellIDs = boundaryCellIDs_[patchI];
        labelList& curBoundaryCellFaceIDs = boundaryCellFaceIDs_[patchI];

        curBoundaryCellIDs.setSize(patchFaces.size());
        curBoundaryCellFaceIDs.setSize(patchFaces.size());

        const labelListList& PointCells = pointCells();

        forAll(patchFaces, faceI)
        {
            bool found = false;

            const face& curFace = patchFaces[faceI];
            const labelList& facePoints = curFace;

            forAll(facePoints, pointI)
            {
                const labelList& facePointCells =
                    PointCells[facePoints[pointI]];

                forAll(facePointCells, cellI)
                {
                    const label curCellIndex = facePointCells[cellI];

                    const faceList& curCellFaces =
                        cellFaces_[curCellIndex];

                    forAll(curCellFaces, cellFaceI)
                    {
                        if (starEqualFace(curFace, curCellFaces[cellFaceI]))
                        {
                            // Found the cell face corresponding to this face
                            found = true;

                            // Set boundary face to the corresponding cell face
                            curBoundaryCellIDs[faceI] = curCellIndex;
                            curBoundaryCellFaceIDs[faceI] = cellFaceI;
                        }
                        if (found) break;
                    }
                    if (found) break;
                }
                if (found) break;
            }
            if (!found)
            {
                FatalErrorIn("starMesh::markBoundaryFaces()")
                    << "Face " << faceI
                    << " does not have neighbour cell."
                    << " Face : " << endl << curFace << endl
                    << "PROSTAR Command: vset,news,vlis";

                forAll(curFace, spI)
                {
                    if (curFace[spI] > -1 && curFace[spI] < starPointID_.size())
                    {
                        Info<< "," << starPointID_[curFace[spI]];
                    }
                    else
                    {
                        Info<< ",???";
                    }
                }

                FatalError
                    << " $ bset,add,vset,all"
                    << abort(FatalError);
            }
        }
    }
}


void starMesh::collectBoundaryFaces()
{
    Info<< "Collecting boundary faces" << endl;
    forAll(boundary_, patchI)
    {
        faceList& patchFaces = boundary_[patchI];

        // set size of patch lists
        const labelList& curBoundaryCellIDs = boundaryCellIDs_[patchI];
        const labelList& curBoundaryCellFaceIDs = boundaryCellFaceIDs_[patchI];

        forAll(curBoundaryCellIDs, faceI)
        {
            patchFaces[faceI] =
                cellFaces_[curBoundaryCellIDs[faceI]]
                    [curBoundaryCellFaceIDs[faceI]];
        }
    }

    Info<< "Finished collecting boundary faces" << endl;
}


// ************************************************************************* //
