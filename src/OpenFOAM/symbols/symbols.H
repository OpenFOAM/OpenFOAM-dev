/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::symbols

Description
    ...

SourceFiles
    ...

\*---------------------------------------------------------------------------*/

#ifndef symbols_H
#define symbols_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace symbols
{

/*---------------------------------------------------------------------------*\
                         Class tokeniser Declaration
\*---------------------------------------------------------------------------*/

class tokeniser
{
    // Private Data

        //- Reference to the stream
        Istream& is_;

        //- List of tokens
        List<token> tokens_;

        //- Index of the first token
        label start_;

        //- Total number of tokens
        label size_;


    // Private Member Functions

        //- Add a token to the end of the stack
        void push(const token&);

        //- Remove a token from the end of the stack
        token pop();

        //- Add a popped token back to the end of the stack
        void unpop(const token&);


public:

    // Constructors

        //- Construct for an input stream
        tokeniser(Istream&);


    // Member Functions

        //- Access the stream
        Istream& stream();

        //- ...
        bool hasToken() const;

        //- ...
        token nextToken();

        //- ...
        void putBack(const token&);

        //- ...
        void splitWord(const word&);

        //- ...
        static bool valid(char c);

        //- ...
        static label priority(const token& t);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- Parse tokens into a dimension set or unit conversion, assuming that the '['
//  has already been read
template<class Type>
Type parseNoBegin
(
    const label lastPrior,
    tokeniser& tis,
    const Type& identity,
    const HashTable<Type>& table
);

//- Parse tokens into a dimension set or unit conversion, assuming that the '['
//  has already been read, and ending before the '[' or ':'. Note that this
//  will leave the stream with a pushed back token.
template<class Type>
Type parseNoBeginOrEnd
(
    Istream& is,
    const Type& identity,
    const HashTable<Type>& table
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace symbols
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "symbolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
