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
    Foam::keyType

Description
    A class for handling keywords in dictionaries.

    A keyType is the keyword of a dictionary.
    It differs from word in that it accepts patterns (regular expressions).

SourceFiles
    keyType.C

\*---------------------------------------------------------------------------*/

#ifndef keyType_H
#define keyType_H

#include "functionName.H"
#include "variable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;
class token;


// Forward declaration of friend functions and operators

class keyType;

Istream& operator>>(Istream&, keyType&);
Ostream& operator<<(Ostream&, const keyType&);


/*---------------------------------------------------------------------------*\
                           Class keyType Declaration
\*---------------------------------------------------------------------------*/

class keyType
:
    public variable
{
    //- Enumeration of the keyword types
    enum type
    {
        UNDEFINED,
        WORD,
        FUNCTIONNAME,
        VARIABLE,
        PATTERN
    };


    // Private Data

        //- The type of this keyword
        type type_;


public:

    // Static Data Members

        //- An empty keyType
        static const keyType null;


    // Constructors

        //- Construct null
        inline keyType();

        //- Copy constructor
        inline keyType(const keyType&);

        //- Construct as word
        inline keyType(const word&);

        //- Construct as functionName
        inline explicit keyType(const functionName&);

        //- Construct as variable
        inline explicit keyType(const variable&);

        //- Construct as pattern
        inline explicit keyType(const string&);

        //- Construct as word from character array.
        inline keyType(const char*);

        //- Construct from token
        explicit keyType(const token&);

        //- Construct from Istream
        explicit keyType(Istream&);


    // Member Functions

        //- Return true if the type has not been defined
        inline bool isUndefined() const;

        //- Return true if the keyword is a functionName
        inline bool isFunctionName() const;

        //- Return true if the keyword is a variable
        inline bool isVariable() const;

        //- Should be treated as a match rather than a literal string
        inline bool isPattern() const;

        //- Smart match as regular expression or as a string
        //  Optionally force a literal match only
        bool match(const std::string&, bool literalMatch=false) const;


    // Member Operators

        // Assignment

            //- Assignment operator
            inline void operator=(const keyType&);

            //- Assign as functionName
            inline void operator=(const functionName&);

            //- Assign as variable
            inline void operator=(const variable&);

            //- Assign as word
            inline void operator=(const word&);

            //- Assign as pattern
            inline void operator=(const string&);

            //- Assign as word
            inline void operator=(const char*);

            //- Disallow assignments where we cannot determine string/word type
            void operator=(const std::string&) = delete;

            //- Assign from token setting the appropriate type
            void operator=(const token&);


    // IOstream Operators

        friend Istream& operator>>(Istream&, keyType&);
        friend Ostream& operator<<(Ostream&, const keyType&);
};


//- Write the keyword to the Ostream with the current level of indentation
Ostream& writeKeyword(Foam::Ostream& os, const keyType& kw);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "keyTypeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
