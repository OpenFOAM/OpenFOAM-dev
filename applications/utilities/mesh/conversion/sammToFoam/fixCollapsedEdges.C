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
    Create intermediate mesh files from SAMM files

\*---------------------------------------------------------------------------*/

#include "sammMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::sammMesh::fixCollapsedEdges()
{
    cellFaces_.setSize(cellShapes_.size());

    forAll(cellShapes_, celli)
    {
        cellFaces_[celli] = cellShapes_[celli].faces();
    }

    // go through the faces and find if there exist faces with duplicate
    // vertices. If so, purge the duplicates and mark the mesh as a polyMesh

    forAll(cellFaces_, celli)
    {
        faceList& curFaces = cellFaces_[celli];

        forAll(curFaces, facei)
        {
            face& vertexLabels = curFaces[facei];

            bool duplicatesFound = false;

            forAll(vertexLabels, vI)
            {
                label curLabel = vertexLabels[vI];

                label nFound = 0;

                forAll(vertexLabels, searchI)
                {
                    if (vertexLabels[searchI] == curLabel)
                    {
                        nFound++;
                    }
                }

                if (nFound > 1)
                {
                    duplicatesFound = true;

                    break;
                }
            }

            if (duplicatesFound)
            {
                // this mesh cannot be described as a shapeMesh
                isShapeMesh_ = false;

                // I am not allowed to reset the shape pointer to unknown
                // here as the shape is still needed to determine which face
                // of the shape is used in potential couple matches.  This
                // will be done in the end using the purgeShapes()
                //

                // create a new face without duplicates and replace original
                face newFace(vertexLabels.size());

                label nNewVertices = 0;

                forAll(vertexLabels, vI)
                {
                    // In order for a face to be a valid entity, duplicate
                    // vertices can only be consecutive (otherwise, the
                    // collapse creates an invalid face). We shall use this
                    // property in the creation of the collapsed face

                    label curLabel = vertexLabels[vI];

                    bool found = false;

                    // search through all vertices from the new face. If the
                    // current label has not been added, add it to the end.
                    for (label searchI = 0; searchI < nNewVertices; searchI++)
                    {
                        if (newFace[searchI] == curLabel)
                        {
                            found = true;

                            break;
                        }
                    }

                    if (!found)
                    {
                        newFace[nNewVertices] = curLabel;
                        nNewVertices++;
                    }
                }

                newFace.setSize(nNewVertices);

                // If the number of non-duplicate labels in the face is less
                // than three, the face has been collapsed in an invalid
                // manner. Error.

                if (nNewVertices < 3)
                {
                    FatalErrorInFunction
                        << "face " << facei << " of cell " << celli
                        << " is collapsed down to a point or edge, which is "
                        << "not permitted" << endl
                        << "original face: " << vertexLabels << endl
                        << "purged face: " << newFace << endl
                        << abort(FatalError);
                }
                else
                {
                    vertexLabels = newFace;
                }
            }
        }
    }
}


// ************************************************************************* //
