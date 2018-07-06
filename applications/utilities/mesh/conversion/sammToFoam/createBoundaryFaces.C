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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::sammMesh::sammEqualFace
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


void Foam::sammMesh::createBoundaryFaces()
{
    forAll(boundary_, patchi)
    {
        faceList& patchFaces = boundary_[patchi];

        const labelListList& PointCells = pointCells();

        forAll(patchFaces, facei)
        {
            bool found = false;

            face& curFace = patchFaces[facei];
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
                        if (sammEqualFace(curCellFaces[cellFacei], curFace))
                        {
                            // Found the cell face corresponding to this face
                            found = true;

                            // Set boundary face to the corresponding cell face
                            // which guarantees it is outward-pointing
                            curFace = curCellFaces[cellFacei];
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
                    << " does not have neighbour cell." << endl
                    << "    face : " << endl << curFace
                    << abort(FatalError);
            }
        }
    }
}


// ************************************************************************* //
