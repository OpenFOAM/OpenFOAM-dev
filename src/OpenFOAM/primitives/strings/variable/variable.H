/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::variable

Description
    A variable is a word with support for additional characters, in particular
    '$' and '/'.

SourceFiles
    variable.C

\*---------------------------------------------------------------------------*/

#ifndef variable_H
#define variable_H

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class variable;
class Istream;
class Ostream;

Istream& operator>>(Istream&, variable&);
Ostream& operator<<(Ostream&, const variable&);


/*---------------------------------------------------------------------------*\
                          Class variable Declaration
\*---------------------------------------------------------------------------*/

class variable
:
    public word
{
    // Private Member Functions

        //- Strip invalid characters from this variable
        inline void stripInvalid();


public:

    // Static Data Members

        static const char* const typeName;
        static int debug;

        //- An empty variable
        static const variable null;


    // Constructors

        //- Construct null
        inline variable();

        //- Copy constructor
        inline variable(const variable&);

        //- Copy constructor of word
        inline explicit variable(const word&);

        //- Copy constructor of string
        inline explicit variable(const string&, const bool doStripInvalid=true);

        //- Copy constructor of std::string
        inline explicit variable
        (
            const std::string&,
            const bool doStripInvalid=true
        );

        //- Copy constructor of character array
        inline explicit variable(const char*, const bool doStripInvalid=true);

        //- Construct from Istream
        //  Words are treated as literals, strings with an auto-test
        explicit variable(Istream&);


    // Member Functions

        //- Is this character valid for a variable
        inline static bool valid(char);


    // Member Operators

        // Assignment

            inline void operator=(const variable&);
            inline void operator=(const word&);
            inline void operator=(const string&);
            inline void operator=(const std::string&);
            inline void operator=(const char*);


    // IOstream Operators

        friend Istream& operator>>(Istream&, variable&);
        friend Ostream& operator<<(Ostream&, const variable&);
};


void writeEntry(Ostream& os, const variable& value);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "variableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
