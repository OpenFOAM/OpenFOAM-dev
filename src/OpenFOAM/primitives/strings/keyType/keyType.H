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

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;


// Forward declaration of friend functions and operators

class keyType;

Istream& operator>>(Istream&, keyType&);
Ostream& operator<<(Ostream&, const keyType&);


/*---------------------------------------------------------------------------*\
                           Class keyType Declaration
\*---------------------------------------------------------------------------*/

class keyType
:
    public word
{
    // Private Data

        //- Is the keyType a pattern (regular expression)
        bool isPattern_;

    // Private Member Functions

        //- Disallow assignments where we cannot determine string/word type
        void operator=(const std::string&);

public:

    // Static Data Members

        //- An empty keyType
        static const keyType null;


    // Constructors

        //- Construct null
        inline keyType();

        //- Copy constructor
        inline keyType(const keyType&);

        //- Copy constructor of word. Not treated as a regular expression
        inline keyType(const word&);

        //- Copy constructor of string. Treat as regular expression.
        inline keyType(const string&);

        //- Copy constructor of character array.
        //  Not treated as a regular expression
        inline keyType(const char*);

        //- Copy constructor of std::string with specified treatment
        inline keyType(const std::string&, const bool isPattern);

        //- Construct from Istream
        //  Treat as regular expression if surrounded by quotation marks.
        keyType(Istream&);


    // Member Functions

        //- Should be treated as a match rather than a literal string
        inline bool isPattern() const;

        //- Smart match as regular expression or as a string
        //  Optionally force a literal match only
        bool match(const std::string&, bool literalMatch=false) const;

    // Member Operators

        // Assignment

            //- Assignment operator
            inline void operator=(const keyType&);

            //- Assign as word, not as non regular expression
            inline void operator=(const word&);

            //- Assign as regular expression
            inline void operator=(const string&);

            //- Assign as word, not as non regular expression
            inline void operator=(const char*);


    // IOstream Operators

        friend Istream& operator>>(Istream&, keyType&);
        friend Ostream& operator<<(Ostream&, const keyType&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "keyTypeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
