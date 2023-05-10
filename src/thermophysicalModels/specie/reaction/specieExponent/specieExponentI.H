/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "specieExponent.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline bool Foam::specieExponent::hasIntegerExponent() const
{
    return integerExponent_ != noIntegerExponent_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::specieExponent::specieExponent()
:
    integerExponent_(noIntegerExponent_),
    scalarExponent_(NaN)
{}


inline Foam::specieExponent::specieExponent(const label integerExponent)
:
    integerExponent_(integerExponent),
    scalarExponent_(integerExponent)
{}


inline Foam::specieExponent::specieExponent(const scalar scalarExponent)
:
    integerExponent_(labelMax),
    scalarExponent_(scalarExponent)
{
    const label integerExponent = floor(scalarExponent);
    if (integerExponent == scalarExponent)
    {
        integerExponent_ = integerExponent;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::specieExponent::operator scalar() const
{
    return scalarExponent_;
}


inline Foam::specieExponent&
Foam::specieExponent::operator=(const label integerExponent)
{
    return *this = specieExponent(integerExponent);
}


inline Foam::specieExponent&
Foam::specieExponent::operator=(const scalar scalarExponent)
{
    return *this = specieExponent(scalarExponent);
}


inline Foam::specieExponent Foam::specieExponent::operator-() const
{
    if (hasIntegerExponent())
    {
        return -integerExponent_;
    }
    else
    {
        return -scalarExponent_;
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::pow(const scalar x, const specieExponent& e)
{
    if (e.hasIntegerExponent())
    {
        return integerPow(x, e.integerExponent_);
    }
    else
    {
        return pow(x, e.scalarExponent_);
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::specieExponent Foam::operator+
(
    const specieExponent& a,
    const specieExponent& b
)
{
    if (a.hasIntegerExponent() && b.hasIntegerExponent())
    {
        return a.integerExponent_ + b.integerExponent_;
    }
    else
    {
        return a.scalarExponent_ + b.scalarExponent_;
    }
}


inline Foam::specieExponent Foam::operator-
(
    const specieExponent& a,
    const specieExponent& b
)
{
    return a + (-b);
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const specieExponent& e)
{
    os  << e.operator scalar();
    return os;
}


// ************************************************************************* //
