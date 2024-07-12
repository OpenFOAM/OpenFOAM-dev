/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::tokenTupleN

Description
    Reads a n-component tuple, storing the result as a list of tokens, using
    bracket counting to disambiguate types. This allows foam-format tables of
    mixed types to be read and columns extracted, provided that the complex
    types are bounded by parentheses. For example, the following table, with
    scalar, vector, scalar, and symmTensor columns could be read as a list of
    tokenTupleN-s:

    \verbatim
        (
            (0 (1 2 3) 4 (5 6 7 8 9 10))
            (1 (2 3 4) 5 (6 7 8 9 10 11))
            (2 (3 4 5) 6 (7 8 9 10 11 12))
        )
    \endverbatim

    Once constructed, this type then provides type-based access to the
    components of the tuple so that they can be extracted into a different
    structure.

SourceFiles
    tokenTupleN.C

\*---------------------------------------------------------------------------*/

#ifndef tokenTupleN_H
#define tokenTupleN_H

#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class tokenTupleN;
Istream& operator>>(Istream&, tokenTupleN&);


/*---------------------------------------------------------------------------*\
                          Class tokenTupleN Declaration
\*---------------------------------------------------------------------------*/

class tokenTupleN
{
    // Private Data

        //- Tokens
        DynamicList<token> tokens_;

        //- Indices of the tokens that start the different pieces of data, plus
        //  an index at the end for convenience
        DynamicList<label> offsets_;


public:

    // Constructors

        //- Default construct
        tokenTupleN();

        //- Construct from stream
        tokenTupleN(Istream& is);


    //- Destructor
    ~tokenTupleN();


    // Member Functions

        //- Get the number of elements in the tuple
        inline label size() const;

        //- Get a type at a given position in the tuple. Stream used for error
        //  reporting only.
        template<class Type>
        inline Type get(const IOstream& s, const label i);


    // IOstream Operators

        //- Read from stream
        friend Istream& operator>>(Istream&, tokenTupleN&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tokenTupleNI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
