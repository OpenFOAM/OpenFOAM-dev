/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include <cctype>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::variable::stripInvalid()
{
    // skip stripping unless debug is active to avoid
    // costly operations
    if (debug && string::stripInvalid<variable>(*this))
    {
        std::cerr
            << "variable::stripInvalid() called for variable "
            << this->c_str() << std::endl;

        if (debug > 1)
        {
            std::cerr
                << "    For debug level (= " << debug
                << ") > 1 this is considered fatal" << std::endl;
            std::abort();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::variable::variable()
:
    word()
{}


inline Foam::variable::variable(const variable& v)
:
    word(v)
{}


inline Foam::variable::variable(const word& w)
:
    word(w)
{}


inline Foam::variable::variable(const string& s, const bool doStripInvalid)
:
    word(s)
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


inline Foam::variable::variable(const std::string& s, const bool doStripInvalid)
:
    word(s)
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


inline Foam::variable::variable(const char* s, const bool doStripInvalid)
:
    word(s)
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::variable::valid(char c)
{
    return
    (
        !isspace(c)
     && c != '"'   // string quote
     && c != '\''  // string quote
     && c != ';'   // end statement
     && c != '{'   // beg subdict
     && c != '}'   // end subdict
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::variable::operator=(const variable& s)
{
    string::operator=(s);
}


inline void Foam::variable::operator=(const word& s)
{
    string::operator=(s);
}


inline void Foam::variable::operator=(const string& s)
{
    string::operator=(s);
    stripInvalid();
}


inline void Foam::variable::operator=(const std::string& s)
{
    string::operator=(s);
    stripInvalid();
}


inline void Foam::variable::operator=(const char* s)
{
    string::operator=(s);
    stripInvalid();
}


// ************************************************************************* //
