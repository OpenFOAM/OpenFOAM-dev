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
    Foam::functionName

Description
    A functionName is a word starting with '#'.

SourceFiles
    functionName.C

\*---------------------------------------------------------------------------*/

#ifndef functionName_H
#define functionName_H

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class functionName;
class Istream;
class Ostream;

Istream& operator>>(Istream&, functionName&);
Ostream& operator<<(Ostream&, const functionName&);


/*---------------------------------------------------------------------------*\
                        Class functionName Declaration
\*---------------------------------------------------------------------------*/

class functionName
:
    public word
{

public:

    // Static Data Members

        static const char* const typeName;
        static int debug;

        //- An empty functionName
        static const functionName null;


    // Constructors

        //- Construct null
        inline functionName();

        //- Copy constructor
        inline functionName(const functionName&);

        //- Copy constructor of word
        inline explicit functionName(const word&);

        //- Copy constructor of string
        inline explicit functionName
        (
            const string&,
            const bool doStripInvalid=true
        );

        //- Copy constructor of std::string
        inline explicit functionName
        (
            const std::string&,
            const bool doStripInvalid=true
        );

        //- Copy constructor of character array
        inline explicit functionName
        (
            const char*,
            const bool doStripInvalid=true
        );

        //- Construct from Istream
        //  Words are treated as literals, strings with an auto-test
        explicit functionName(Istream&);


    // Member Functions

        //- Is this character valid for a functionName
        inline static bool valid(char);


    // Member Operators

        // Assignment

            inline void operator=(const functionName&);
            inline void operator=(const word&);
            inline void operator=(const string&);
            inline void operator=(const std::string&);
            inline void operator=(const char*);


    // IOstream Operators

        friend Istream& operator>>(Istream&, functionName&);
        friend Ostream& operator<<(Ostream&, const functionName&);
};


void writeEntry(Ostream& os, const functionName& value);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "functionNameI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
