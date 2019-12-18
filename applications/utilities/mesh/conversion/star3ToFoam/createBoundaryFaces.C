/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
// PROSTAR boundary format idiosyncrasies
bool Foam::starMesh::starEqualFace
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


void Foam::starMesh::markBoundaryFaces()
{
    // set size of mark lists for the boundary
    boundaryCellIDs_.setSize(boundary_.size());
    boundaryCellFaceIDs_.setSize(boundary_.size());

    forAll(boundary_, patchi)
    {
        const faceList& patchFaces = boundary_[patchi];

        // set size of patch lists
        labelList& curBoundaryCellIDs = boundaryCellIDs_[patchi];
        labelList& curBoundaryCellFaceIDs = boundaryCellFaceIDs_[patchi];

        curBoundaryCellIDs.setSize(patchFaces.size());
        curBoundaryCellFaceIDs.setSize(patchFaces.size());

        const labelListList& PointCells = pointCells();

        forAll(patchFaces, facei)
        {
            bool found = false;

            const face& curFace = patchFaces[facei];
            const labelList& facePoints = curFace;

            forAll(facePoints, pointi)
            {
                const labelList& facePointCells =
                    PointCells[facePoints[pointi]];

                forAll(facePointCells, celli)
                {
                    const label curCellIndex = facePointCells[celli];

                    const faceList& curCellFaces =
                        cellFaces_[curCellIndex];

                    forAll(curCellFaces, cellFacei)
                    {
                        if (starEqualFace(curFace, curCellFaces[cellFacei]))
                        {
                            // Found the cell face corresponding to this face
                            found = true;

                            // Set boundary face to the corresponding cell face
                            curBoundaryCellIDs[facei] = curCellIndex;
                            curBoundaryCellFaceIDs[facei] = cellFacei;
                        }
                        if (found) break;
                    }
                    if (found) break;
                }
                if (found) break;
            }
            if (!found)
            {
                FatalErrorInFunction
                    << "Face " << facei
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


void Foam::starMesh::collectBoundaryFaces()
{
    Info<< "Collecting boundary faces" << endl;
    forAll(boundary_, patchi)
    {
        faceList& patchFaces = boundary_[patchi];

        // set size of patch lists
        const labelList& curBoundaryCellIDs = boundaryCellIDs_[patchi];
        const labelList& curBoundaryCellFaceIDs = boundaryCellFaceIDs_[patchi];

        forAll(curBoundaryCellIDs, facei)
        {
            patchFaces[facei] =
                cellFaces_[curBoundaryCellIDs[facei]]
                    [curBoundaryCellFaceIDs[facei]];
        }
    }

    Info<< "Finished collecting boundary faces" << endl;
}


// ************************************************************************* //
