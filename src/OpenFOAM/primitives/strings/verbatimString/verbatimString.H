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
    Foam::verbatimString

Description
    A class for handling verbatimStrings, derived from string.

    A verbatimString is a verbatim string of characters.

SourceFiles
    verbatimString.C
    verbatimStringIO.C

\*---------------------------------------------------------------------------*/

#ifndef verbatimString_H
#define verbatimString_H

#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class verbatimString;
inline verbatimString operator&(const verbatimString&, const verbatimString&);
Istream& operator>>(Istream&, verbatimString&);
Ostream& operator<<(Ostream&, const verbatimString&);


/*---------------------------------------------------------------------------*\
                       Class verbatimString Declaration
\*---------------------------------------------------------------------------*/

class verbatimString
:
    public string
{

public:

    // Static Data Members

        static const char* const typeName;
        static int debug;

        //- An empty verbatimString
        static const verbatimString null;


    // Constructors

        //- Construct null
        inline verbatimString();

        //- Copy constructor
        inline verbatimString(const verbatimString&);

        //- Copy constructor of character array
        inline verbatimString(const char*);

        //- Copy constructor of string
        inline verbatimString(const string&);

        //- Copy constructor of std::string
        inline verbatimString(const std::string&);

        //- Construct from Istream
        verbatimString(Istream&);


    // Member Operators

        // Assignment

            inline void operator=(const verbatimString&);
            inline void operator=(const string&);
            inline void operator=(const std::string&);
            inline void operator=(const char*);


    // IOstream Operators

        friend Istream& operator>>(Istream&, verbatimString&);
        friend Ostream& operator<<(Ostream&, const verbatimString&);
};


void writeEntry(Ostream& os, const verbatimString& value);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "verbatimStringI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
