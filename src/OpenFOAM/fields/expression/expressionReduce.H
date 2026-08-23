/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

Classes
    Foam::expression::Reduce
    Foam::expression::ReduceIterable

Description
    Class which provides the ability to reduce an expression to a single value
    using a binary operation (probably either min, max or sum). This must be
    specialised for expressionable containers.

    In the case of iterable containers that specialisation can inherit a
    generic iterable implementation of reduction from the ReduceIterable class
    (also defined in this file).

    Non-iterable containers must implement their own reduction methods. This is
    likely to be straightforward as non-iterable containers are typically just
    wrappings or structured collections of iterable containers, so the
    non-iterable implementation can just forward to the iterable implementation
    and combine and return the result.

\*---------------------------------------------------------------------------*/

#ifndef expressionReduce_H
#define expressionReduce_H

#include "expressionResult.H"
#include <numeric>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

/*---------------------------------------------------------------------------*\
                          Class Reduce Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct Reduce
{
    static const bool globalReduce = false;

    template<class Op, class Value>
    static const Value& value(const Value& x)
    {
        return x;
    }
};

/*---------------------------------------------------------------------------*\
                      Class ReduceIterable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct ReduceIterable;

template<template<class ...> class Container, class Type, class ... Args>
struct ReduceIterable<Container<Type, Args ...>>
{
    static const bool globalReduce = false;

    template<class Op, class Expression>
    static typename ResultIterable<Container<Type, Args ...>>::
    template elementType<Expression> value(const Expression& e)
    {
        typedef typename ResultIterable<Container<Type, Args ...>>::template
            operateType<Expression> operateType;

        typedef typename ResultIterable<Container<Type, Args ...>>::template
            elementType<Expression> elementType;

        return
            std::accumulate
            (
                begin(e),
                end(e),
                Op::template identity<elementType>(),
                [](const elementType& x, const operateType& y)
                {
                    return Op::value(x, Reduce<Type>::template value<Op>(y));
                }
            );
    }
};

/*---------------------------------------------------------------------------*\
                    Helper functions for primitive reductions
\*---------------------------------------------------------------------------*/

template<class Type>
inline Type getPrimitiveValue(const Type& s)
{
    return s;
}

template<class Type>
inline void setPrimitiveValue(Type& s, const Type& ps)
{
    s = ps;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Expression>
using ReduceType = typename expression::Reduce
<
    typename expression::AllContainers<Expression>::type
>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
