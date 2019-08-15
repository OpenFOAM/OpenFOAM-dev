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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::keyType::keyType()
:
    variable(),
    isPattern_(false)
{}


inline Foam::keyType::keyType(const keyType& k)
:
    variable(k),
    isPattern_(k.isPattern())
{}


inline Foam::keyType::keyType(const word& w)
:
    variable(w),
    isPattern_(false)
{}


inline Foam::keyType::keyType(const variable& v)
:
    variable(v),
    isPattern_(false)
{}


inline Foam::keyType::keyType(const string& s)
:
    variable(s, false),
    isPattern_(true)
{}


inline Foam::keyType::keyType(const char* s)
:
    variable(s, false),
    isPattern_(false)
{}


inline Foam::keyType::keyType(const std::string& s, const bool isPattern)
:
    variable(s, false),
    isPattern_(isPattern)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::keyType::isPattern() const
{
    return isPattern_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::keyType::operator=(const keyType& k)
{
    variable::operator=(k);
    isPattern_ = k.isPattern_;
}


inline void Foam::keyType::operator=(const variable& v)
{
    variable::operator=(v);
    isPattern_ = false;
}


inline void Foam::keyType::operator=(const word& w)
{
    variable::operator=(w);
    isPattern_ = false;
}


inline void Foam::keyType::operator=(const string& s)
{
    // Bypass checking
    string::operator=(s);
    isPattern_ = true;
}


inline void Foam::keyType::operator=(const char* s)
{
    // Bypass checking
    string::operator=(s);
    isPattern_ = false;
}


// ************************************************************************* //
