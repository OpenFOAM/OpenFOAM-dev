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

\*---------------------------------------------------------------------------*/

#include "zeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::scalar Foam::zeroField::operator[](const label) const
{
    return scalar(0);
}


inline Foam::zeroField Foam::zeroField::primitiveField() const
{
    return zeroField();
}


inline Foam::zeroField Foam::zeroField::operator()() const
{
    return zeroField();
}


inline Foam::zeroField Foam::zeroField::operator-() const
{
    return zeroField();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const zeroField operator+(const zeroField&, const zeroField&)
{
    return zeroField();
}

template<class Type>
inline const Type& operator+(const Type& t, const zeroField&)
{
    return t;
}

template<class Type>
inline const Type& operator+(const zeroField&, const Type& t)
{
    return t;
}

inline const zeroField operator-(const zeroField&, const zeroField&)
{
    return zeroField();
}

template<class Type>
inline const Type& operator-(const Type& t, const zeroField&)
{
    return t;
}

template<class Type>
inline Type operator-(const zeroField&, const Type& t)
{
    return -t;
}

template<class Type>
inline zeroField operator*(const Type& t, const zeroField&)
{
    return zeroField();
}

template<class Type>
inline zeroField operator*(const zeroField&, const Type& t)
{
    return zeroField();
}

template<class Type>
inline zeroField operator/(const zeroField&, const Type& t)
{
    return zeroField();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
