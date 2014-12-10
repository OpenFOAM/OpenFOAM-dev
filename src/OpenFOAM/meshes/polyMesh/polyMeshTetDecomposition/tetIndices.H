/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "tetrahedron.H"
#include "triPointRef.H"
#include "polyMesh.H"
#include "triFace.H"
#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class tetIndices Declaration
\*---------------------------------------------------------------------------*/

class tetIndices
{
    // Private data

    //- cell that this is a decomposed tet of
    label cellI_;

    //- face that holds this decomposed tet
    label faceI_;

    //- base point on the face
    label faceBasePtI_;

    //- point on the face such that the right-hand circulation
    //      {faceBasePtI_, facePtAI_, facePtBI_}
    //  forms a triangle that points out of the tet
    label facePtAI_;

    //- point on the face such that the right-hand circulation
    //      {faceBasePtI_, facePtAI_, facePtBI_}
    //  forms a triangle that points out of the tet
    label facePtBI_;

    //- point on the face, *relative to the base point*, which
    //  characterises this tet on the face.
    label tetPtI_;


public:

    // Constructors

        //- Construct null
        tetIndices();

        //- Construct from components
        tetIndices
        (
            label cellI,
            label faceI,
            label faceBasePtI,
            label facePtAI,
            label facePtBI,
            label tetPtI
        );

        //- Construct from cell, face, pt description of tet
        tetIndices
        (
            label cellI,
            label faceI,
            label tetPtI,
            const polyMesh& mesh
        );


    //- Destructor
    ~tetIndices();


    // Member Functions

        // Access

            //- Return the cell
            inline label cell() const;

            //- Return the face
            inline label face() const;

            //- Return the face base point
            inline label faceBasePt() const;

            //- Return face point A
            inline label facePtA() const;

            //- Return face point B
            inline label facePtB() const;

            //- Return the characterising tetPtI
            inline label tetPt() const;

            //- Return the geometry corresponding to this tet from the
            //  supplied mesh
            inline tetPointRef tet(const polyMesh& mesh) const;

            //- Return the geometry corresponding to this tet from the
            //  supplied mesh using the old positions
            inline tetPointRef oldTet(const polyMesh& mesh) const;

            //- Return the geometry corresponding to the tri on the
            //  mesh face for this tet from the supplied mesh. Normal of
            //  the tri points out of the cell.
            inline triPointRef faceTri(const polyMesh& mesh) const;

            //- Return the point indices corresponding to the tri on the mesh
            //  face for this tet from the supplied mesh. Normal of
            //  the tri points out of the cell.
            inline triFace faceTriIs(const polyMesh& mesh) const;

            //- Return the geometry corresponding to the tri on the
            //  mesh face for this tet from the supplied mesh using
            //  the old position
            inline triPointRef oldFaceTri(const polyMesh& mesh) const;


        // Edit

            //- Return non-const access to the cell
            inline label& cell();

            //- Return non-const access to the face
            inline label& face();

            //- Return non-const access to the face base point
            inline label& faceBasePt();

            //- Return non-const access to face point A
            inline label& facePtA();

            //- Return non-const access to face point B
            inline label& facePtB();

            //- Return non-const access to the characterising tetPtI
            inline label& tetPt();


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
