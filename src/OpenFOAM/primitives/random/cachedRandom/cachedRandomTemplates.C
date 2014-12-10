/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "cachedRandom.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::cachedRandom::sample01()
{
    Type value;
    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        value.component(cmpt) = scalar01();
    }

    return value;
}


template<class Type>
Type Foam::cachedRandom::position(const Type& start, const Type& end)
{
    Type value(start);
    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        value.component(cmpt) +=
            scalar01()*(end.component(cmpt) - start.component(cmpt));
    }

    return value;
}


template<class Type>
void Foam::cachedRandom::randomise01(Type& value)
{
    value = sample01<Type>();
}


template<class Type>
Type Foam::cachedRandom::globalSample01()
{
    Type value = -GREAT*pTraits<Type>::one;

    if (Pstream::master())
    {
        value = sample01<Type>();
    }

    reduce(value, maxOp<Type>());

    return value;
}


template<class Type>
Type Foam::cachedRandom::globalPosition(const Type& start, const Type& end)
{
    Type value = -GREAT*pTraits<Type>::one;

    if (Pstream::master())
    {
        value = position<Type>(start, end);
    }

    reduce(value, maxOp<Type>());

    return value;
}


template<class Type>
void Foam::cachedRandom::globalRandomise01(Type& value)
{
    value = -GREAT*pTraits<Type>::one;

    if (Pstream::master())
    {
        value = sample01<Type>();
    }

    reduce(value, maxOp<Type>());
}


// ************************************************************************* //
