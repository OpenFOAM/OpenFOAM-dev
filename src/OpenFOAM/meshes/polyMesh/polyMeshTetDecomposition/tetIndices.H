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

Class
    Foam::tetIndices

Description
    Storage and named access for the indices of a tet which is part of
    the decomposition of a cell.

    Tets are designated by
    - cell (of course)
    - face on cell
    - three points on face (faceBasePt, facePtA, facePtB)
      When constructing from a mesh and index in the face (tetPtI):
        - faceBasePt is the mesh.tetBasePtIs() base point
        - facePtA is tetPtI away from faceBasePt
        - facePtB is next one after/before facePtA
        e.g.:

            +---+
            |2 /|
            | / |
            |/ 1|  <- tetPt (so 1 for first triangle, 2 for second)
            +---+
            ^
           faceBasePt

SourceFiles
    tetIndicesI.H
    tetIndices.C

\*---------------------------------------------------------------------------*/

#ifndef tetIndices_H
#define tetIndices_H

#include "label.H"
#include "tetPointRef.H"
#include "triPointRef.H"
#include "polyMesh.H"
#include "triFace.H"
#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class tetIndices;

Istream& operator>>(Istream&, tetIndices&);
Ostream& operator<<(Ostream&, const tetIndices&);


/*---------------------------------------------------------------------------*\
                         Class tetIndices Declaration
\*---------------------------------------------------------------------------*/

class tetIndices
{
    // Private Data

        //- Cell that this is a decomposed tet of
        label celli_;

        //- Face that holds this decomposed tet
        label facei_;

        //- Point on the face, *relative to the base point*, which
        //  characterises this tet on the face.
        label tetPti_;


public:

    // Constructors

        //- Construct null
        tetIndices();

        //- Construct from components
        tetIndices(label celli, label facei, label tetPtI);


    //- Destructor
    ~tetIndices();


    // Member Functions

        // Access

            //- Return the cell
            inline label cell() const;

            //- Return non-const access to the cell
            inline label& cell();

            //- Return the face
            inline label face() const;

            //- Return non-const access to the face
            inline label& face();

            //- Return the characterising tetPtI
            inline label tetPt() const;

            //- Return non-const access to the characterising tetPtI
            inline label& tetPt();

            //- Return the indices corresponding to the tri on the face for
            //  this tet. The normal of the tri points out of the cell.
            inline triFace faceTriIs(const polyMesh& mesh) const;

            //- Return the geometry corresponding to this tet
            inline tetPointRef tet(const polyMesh& mesh) const;

            //- Return the geometry corresponding to the tri on the face for
            //  this tet. The normal of the tri points out of the cell.
            inline triPointRef faceTri(const polyMesh& mesh) const;

            //- Return the geometry corresponding to the tri on the face for
            //  this tet using the old positions.
            inline triPointRef oldFaceTri(const polyMesh& mesh) const;


    // Member Operators

        inline bool operator==(const tetIndices&) const;
        inline bool operator!=(const tetIndices&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, tetIndices&);
        friend Ostream& operator<<(Ostream&, const tetIndices&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetIndicesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
