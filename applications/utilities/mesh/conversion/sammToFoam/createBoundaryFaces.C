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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialist version of face comparison to deal with
// PROSTAR boundary format idiosyncracies
bool sammMesh::sammEqualFace
(
    const face& boundaryFace,
    const face& cellFace
) const
{
    // A PROSTAR boundary face is defined by 4 vertices irrespective
    // of its topology.
    // In order to deal with all possibilities, two faces will be
    // considered equal if three of the vertices are the same.
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
        return true;
    }
    else
    {
        return false;
    }
}


void sammMesh::createBoundaryFaces()
{
    forAll(boundary_, patchI)
    {
        faceList& patchFaces = boundary_[patchI];

        const labelListList& PointCells = pointCells();

        forAll(patchFaces, faceI)
        {
            bool found = false;

            face& curFace = patchFaces[faceI];
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
                        if (sammEqualFace(curCellFaces[cellFaceI], curFace))
                        {
                            // Found the cell face corresponding to this face
                            found = true;

                            // Set boundary face to the corresponding cell face
                            // which guarantees it is outward-pointing
                            curFace = curCellFaces[cellFaceI];
                        }
                        if (found) break;
                    }
                    if (found) break;
                }
                if (found) break;
            }
            if (!found)
            {
                FatalErrorIn("sammMesh::createBoundaryFaces()")
                    << "Face " << faceI
                    << " does not have neighbour cell." << endl
                    << "    face : " << endl << curFace
                    << abort(FatalError);
            }
        }
    }
}


// ************************************************************************* //
