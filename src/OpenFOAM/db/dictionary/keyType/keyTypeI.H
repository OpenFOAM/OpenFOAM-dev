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
    type_(UNDEFINED)
{}


inline Foam::keyType::keyType(const keyType& k)
:
    variable(k),
    type_(k.type_)
{}


inline Foam::keyType::keyType(const word& w)
:
    variable(w),
    type_(WORD)
{}


inline Foam::keyType::keyType(const functionName& fn)
:
    variable(fn),
    type_(FUNCTIONNAME)
{}


inline Foam::keyType::keyType(const variable& v)
:
    variable(v),
    type_(VARIABLE)
{}


inline Foam::keyType::keyType(const string& s)
:
    variable(s, false),
    type_(PATTERN)
{}


inline Foam::keyType::keyType(const char* s)
:
    variable(s, true),
    type_(WORD)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::keyType::isUndefined() const
{
    return type_ == UNDEFINED;
}


inline bool Foam::keyType::isFunctionName() const
{
    return type_ == FUNCTIONNAME;
}


inline bool Foam::keyType::isVariable() const
{
    return type_ == VARIABLE;
}


inline bool Foam::keyType::isPattern() const
{
    return type_ == PATTERN;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::keyType::operator=(const keyType& k)
{
    variable::operator=(k);
    type_ = k.type_;
}


inline void Foam::keyType::operator=(const word& w)
{
    variable::operator=(w);
    type_ = WORD;
}


inline void Foam::keyType::operator=(const functionName& fn)
{
    word::operator=(fn);
    type_ = FUNCTIONNAME;
}


inline void Foam::keyType::operator=(const variable& v)
{
    variable::operator=(v);
    type_ = VARIABLE;
}


inline void Foam::keyType::operator=(const string& s)
{
    string::operator=(s);
    type_ = PATTERN;
}


inline void Foam::keyType::operator=(const char* s)
{
    word::operator=(s);
    type_ = WORD;
}


// ************************************************************************* //
