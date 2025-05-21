/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::TypeSet

Description
    Template meta-programming for operations involving sets of types

\*---------------------------------------------------------------------------*/

#ifndef TypeSet_H
#define TypeSet_H

#include "autoPtr.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class TypeSet Declaration
\*---------------------------------------------------------------------------*/

template<class ... Types>
struct TypeSet;

template<class Type, class ... Types>
struct TypeSet<Type, Types ...>
{
    template<class OtherType>
    static bool isAny(const OtherType& t)
    {
        return isA<Type>(t) || TypeSet<Types ...>::isAny(t);
    }

    template<class OtherType>
    static bool isAll(const OtherType& t)
    {
        return isA<Type>(t) && TypeSet<Types ...>::isAll(t);
    }

    //- Clone into the first matching type, or return an empty pointer if a
    //  matching type is not in the set
    template<class OtherType>
    static autoPtr<OtherType> clone(const OtherType& t)
    {
        return
            isA<Type>(t)
          ? autoPtr<OtherType>(new Type(refCast<const Type>(t)))
          : TypeSet<Types ...>::clone(t);
    }
};

template<>
struct TypeSet<>
{
    template<class OtherType>
    static bool isAny(const OtherType&)
    {
        return false;
    }

    template<class OtherType>
    static bool isAll(const OtherType&)
    {
        return true;
    }

    template<class OtherType>
    static autoPtr<OtherType> clone(const OtherType&)
    {
        return autoPtr<OtherType>();
    }
};


/*---------------------------------------------------------------------------*\
                        Class TypeSetContains Declaration
\*---------------------------------------------------------------------------*/

template<class Set, class OtherType>
struct TypeSetContains;

template<class Type, class ... Types, class OtherType>
struct TypeSetContains<TypeSet<Type, Types ...>, OtherType>
:
    public TypeSetContains<TypeSet<Types ...>, OtherType>
{};

template<class Type, class ... Types>
struct TypeSetContains<TypeSet<Type, Types ...>, Type>
:
    public std::true_type
{};

template<class OtherType>
struct TypeSetContains<TypeSet<>, OtherType>
:
    public std::false_type
{};


/*---------------------------------------------------------------------------*\
                      Class TypeSetConcatenate Declaration
\*---------------------------------------------------------------------------*/

template<class SetA, class SetB>
struct TypeSetConcatenate;

template<class ... TypesA, class ... TypesB>
struct TypeSetConcatenate<TypeSet<TypesA ...>, TypeSet<TypesB ...>>
{
    typedef TypeSet<TypesA ..., TypesB ...> type;
};


/*---------------------------------------------------------------------------*\
                        Class TypeSetRemove Declaration
\*---------------------------------------------------------------------------*/

template<class Set, class RemoveSet>
struct TypeSetRemove;

template<class Type, class ... Types, class RemoveSet>
struct TypeSetRemove<TypeSet<Type, Types ...>, RemoveSet>
{
    typedef typename
        TypeSetRemove<TypeSet<Types ...>, RemoveSet>::type
        subType;

    typedef typename
        std::conditional
        <
            TypeSetContains<RemoveSet, Type>::value,
            subType,
            typename TypeSetConcatenate<TypeSet<Type>, subType>::type
        >::type type;
};

template<class RemoveSet>
struct TypeSetRemove<TypeSet<>, RemoveSet>
{
    typedef TypeSet<> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
