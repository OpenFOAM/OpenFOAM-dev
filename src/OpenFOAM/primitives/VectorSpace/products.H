/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Traits classes for inner and outer products of primitives.

\*---------------------------------------------------------------------------*/

#ifndef products_H
#define products_H

#include "pTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt, int rank>
class typeOfRank
{};


template<class Cmpt, int rank>
class symmTypeOfRank
{};


template<class arg1, class arg2>
class typeOfSum
{
public:

    typedef arg1 type;
};


template<class arg1, class arg2>
class outerProduct
{
public:

    typedef typename typeOfRank
    <
        typename pTraits<arg1>::cmptType,
        int(pTraits<arg1>::rank) + int(pTraits<arg2>::rank)
    >::type type;
};


template<class arg1, class arg2>
class crossProduct
{
public:

    typedef typename typeOfRank
    <
        typename pTraits<arg2>::cmptType,
        int(pTraits<arg1>::rank) + int(pTraits<arg2>::rank) - 1
    >::type type;
};

template<class arg1, class arg2>
class innerProduct
{
public:

    typedef typename typeOfRank
    <
        typename pTraits<arg1>::cmptType,
        int(pTraits<arg1>::rank) + int(pTraits<arg2>::rank) - 2
    >::type type;
};

template<class arg1, class arg2>
class scalarProduct
{
public:

    typedef typename pTraits<arg1>::cmptType type;
};


template<class arg1, int arg2>
class powProduct
{
public:

    typedef typename symmTypeOfRank
    <
        typename pTraits<arg1>::cmptType,
        arg2*int(pTraits<arg1>::rank)
    >::type type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
