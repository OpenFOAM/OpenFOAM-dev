/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "Product2.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, Foam::direction rank>
inline Type Foam::Function2s::ProductFunction1s<Type, rank>::value
(
    const scalar x,
    const scalar y
) const
{
    const ProductFilter<Type> filter;
    static const direction otherRank = pTraits<Type>::rank - rank;
    const ProductFunction1s<Type, otherRank>& other = *this;

    if (fs.first().valid())
    {
        return filter(fs.first()->value(x)*other.fs.second()->value(y));
    }
    else if (fs.second().valid())
    {
        return filter(other.fs.first()->value(x)*fs.second()->value(y));
    }
    else
    {
        return ProductFunction1s<Type, rank - 1>::value(x, y);
    }
}


template<class Type>
inline Type Foam::Function2s::ProductFunction1s<Type, 0>::value
(
    const scalar x,
    const scalar y
) const
{
    const ProductFilter<Type> filter;
    return filter(fs.first()->value(x)*fs.second()->value(y));
}


template<class Type>
inline Type Foam::Function2s::Product<Type>::value
(
    const scalar x,
    const scalar y
) const
{
    return fs_.value(x, y);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

inline Foam::word Foam::Function2s::valueName(const direction argument)
{
    return "value" + Foam::name(argument + 1);
}


inline Foam::word Foam::Function2s::valueName
(
    const direction argument,
    const word& typeName
)
{
    return
        typeName.empty()
      ? valueName(argument)
      : valueName(argument) + "<" + typeName + ">";
}


template<class Type>
inline Foam::word Foam::Function2s::valueName(const direction argument)
{
    return valueName(argument, pTraits<Type>::typeName);
}


inline Foam::word Foam::Function2s::valueName
(
    const direction argument,
    const Tuple2<word, label>& typeAndRank
)
{
    return valueName(argument, typeAndRank.first());
}


// ************************************************************************* //
