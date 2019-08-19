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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::functionName::functionName()
:
    word()
{}


inline Foam::functionName::functionName(const functionName& v)
:
    word(v)
{}


inline Foam::functionName::functionName(const word& w)
:
    word(w)
{}


inline Foam::functionName::functionName
(
    const string& s,
    const bool doStripInvalid
)
:
    word(s, doStripInvalid)
{}


inline Foam::functionName::functionName
(
    const std::string& s,
    const bool doStripInvalid
)
:
    word(s, doStripInvalid)
{}


inline Foam::functionName::functionName
(
    const char* s,
    const bool doStripInvalid
)
:
    word(s, doStripInvalid)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::functionName::operator=(const functionName& s)
{
    word::operator=(s);
}


inline void Foam::functionName::operator=(const word& s)
{
    word::operator=(s);
}


inline void Foam::functionName::operator=(const string& s)
{
    word::operator=(s);
}


inline void Foam::functionName::operator=(const std::string& s)
{
    word::operator=(s);
}


inline void Foam::functionName::operator=(const char* s)
{
    word::operator=(s);
}


// ************************************************************************* //
