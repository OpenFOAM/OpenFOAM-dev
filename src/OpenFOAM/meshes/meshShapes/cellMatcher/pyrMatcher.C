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

\*---------------------------------------------------------------------------*/

#include "pyrMatcher.H"
#include "cellMatcher.H"
#include "primitiveMesh.H"
#include "primitiveMesh.H"
#include "cellModeller.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::pyrMatcher::vertPerCell = 5;
const Foam::label Foam::pyrMatcher::facePerCell = 5;
const Foam::label Foam::pyrMatcher::maxVertPerFace = 4;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pyrMatcher::pyrMatcher()
:
    cellMatcher
    (
        vertPerCell,
        facePerCell,
        maxVertPerFace,
        "pyr"
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pyrMatcher::~pyrMatcher()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pyrMatcher::matchShape
(
    const bool checkOnly,
    const faceList& faces,
    const labelList& owner,
    const label celli,
    const labelList& myFaces
)
{
    if (!faceSizeMatch(faces, myFaces))
    {
        return false;
    }

    // Is pyr for sure since no other shape with 1 quad, 4 triangles
    if (checkOnly)
    {
        return true;
    }

    // Calculate localFaces_ and mapping pointMap_, faceMap_
    label numVert = calcLocalFaces(faces, myFaces);

    if (numVert != vertPerCell)
    {
        return false;
    }

    // Set up 'edge' to face mapping.
    calcEdgeAddressing(numVert);

    // Set up point on face to index-in-face mapping
    calcPointFaceIndex();

    // Storage for maps -vertex to mesh and -face to mesh
    vertLabels_.setSize(vertPerCell);
    faceLabels_.setSize(facePerCell);

    //
    // Start from quad face (face0)
    //

    label face0I = -1;
    forAll(faceSize_, facei)
    {
        if (faceSize_[facei] == 4)
        {
            face0I = facei;
            break;
        }
    }
    const face& face0 = localFaces_[face0I];
    label face0vert0 = 0;


    //
    // Try to follow prespecified path on faces of cell,
    // starting at face0vert0
    //

    vertLabels_[0] = pointMap_[face0[face0vert0]];
    faceLabels_[0] = faceMap_[face0I];

    // Walk face 0 from vertex 0 to 1
    label face0vert1 =
        nextVert
        (
            face0vert0,
            faceSize_[face0I],
            !(owner[faceMap_[face0I]] == celli)
        );
    vertLabels_[1] = pointMap_[face0[face0vert1]];

    // Walk face 0 from vertex 1 to 2
    label face0vert2 =
        nextVert
        (
            face0vert1,
            faceSize_[face0I],
            !(owner[faceMap_[face0I]] == celli)
        );
    vertLabels_[2] = pointMap_[face0[face0vert2]];

    // Walk face 0 from vertex 2 to 3
    label face0vert3 =
        nextVert
        (
            face0vert2,
            faceSize_[face0I],
            !(owner[faceMap_[face0I]] == celli)
        );
    vertLabels_[3] = pointMap_[face0[face0vert3]];

    // Jump edge from face0 to face1
    label face1I =
        otherFace
        (
            numVert,
            face0[face0vert3],
            face0[face0vert0],
            face0I
        );
    faceLabels_[1] = faceMap_[face1I];

    // Jump edge from face0 to face2
    label face2I =
        otherFace
        (
            numVert,
            face0[face0vert2],
            face0[face0vert3],
            face0I
        );
    faceLabels_[2] = faceMap_[face2I];

    // Jump edge from face0 to face3
    label face3I =
        otherFace
        (
            numVert,
            face0[face0vert1],
            face0[face0vert2],
            face0I
        );
    faceLabels_[3] = faceMap_[face3I];

    // Jump edge from face0 to face4
    label face4I =
        otherFace
        (
            numVert,
            face0[face0vert0],
            face0[face0vert1],
            face0I
        );
    faceLabels_[4] = faceMap_[face4I];

    const face& face4 = localFaces_[face4I];

    // Get index of vert0 in face 4
    label face4vert0 = pointFaceIndex_[face0[face0vert0]][face4I];

    // Walk face 4 from vertex 0 to 4
    label face4vert4 =
        nextVert
        (
            face4vert0,
            faceSize_[face4I],
            !(owner[faceMap_[face4I]] == celli)
        );
    vertLabels_[4] = pointMap_[face4[face4vert4]];

    return true;
}


Foam::label Foam::pyrMatcher::faceHashValue() const
{
    return 4*3+4;
}


bool Foam::pyrMatcher::faceSizeMatch
(
    const faceList& faces,
    const labelList& myFaces
) const
{
    if (myFaces.size() != 5)
    {
        return false;
    }

    label nTris = 0;
    label nQuads = 0;

    forAll(myFaces, myFacei)
    {
        label size = faces[myFaces[myFacei]].size();

        if (size == 3)
        {
            nTris++;
        }
        else if (size == 4)
        {
            nQuads++;
        }
        else
        {
            return false;
        }
    }

    if ((nTris == 4) && (nQuads == 1))
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::pyrMatcher::isA(const primitiveMesh& mesh, const label celli)
{
    return matchShape
    (
        true,
        mesh.faces(),
        mesh.faceOwner(),
        celli,
        mesh.cells()[celli]
    );
}


bool Foam::pyrMatcher::isA(const faceList& faces)
{
    // Do as if mesh with one cell only
    return matchShape
    (
        true,
        faces,                      // all faces in mesh
        labelList(faces.size(), 0), // cell 0 is owner of all faces
        0,                          // cell label
        identity(faces.size())      // faces of cell 0
    );
}


bool Foam::pyrMatcher::matches
(
    const primitiveMesh& mesh,
    const label celli,
    cellShape& shape
)
{
    if
    (
        matchShape
        (
            false,
            mesh.faces(),
            mesh.faceOwner(),
            celli,
            mesh.cells()[celli]
        )
    )
    {
        shape = cellShape(model(), vertLabels());

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
