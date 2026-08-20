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
    Foam::expression::Result
    Foam::expression::ResultIterable

Description
    Class which provides the types returned by expression evaluations. This
    must be specialised for expressionable containers.

    In the case of iterable containers that specialisation can inherit generic
    definitions from the ResultIterable class (also defined in this file).

    Non-iterable containers must implement their own result classes. This is
    likely to be straightforward as non-iterable containers are typically just
    wrappings or structured collections of iterable containers, so the
    non-iterable implementation can just wrap the types provided by the
    iterable classes.

\*---------------------------------------------------------------------------*/

#ifndef expressionResult_H
#define expressionResult_H

#include "expressionIterator.H"
#include "expressionContainers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

/*---------------------------------------------------------------------------*\
                          Class Result Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct Result
{
    template<class Value>
    using elementType = Value;

    template<class Value>
    using type = Value;
};

/*---------------------------------------------------------------------------*\
                      Class ResultIterable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct ResultIterable;

template<template<class ...> class Container, class Type, class ... Args>
struct ResultIterable<Container<Type, Args ...>>
{
    template<class Expression>
    using operateType =
        std::decay_t
        <
            typename Operate
            <
                std::decay_t<Expression>,
                Arguments
                <
                    typename Operate
                    <
                        std::decay_t<Expression>,
                        Arguments<Bound<false>>
                    >::type
                >
            >::type
        >;

    template<class Expression>
    using elementType =
        typename Result<Type>::template elementType<operateType<Expression>>;

    template<class Expression>
    using type =
        Container
        <
            typename Result<Type>::template type<operateType<Expression>>,
            Args ...
        >;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Expression>
using ElementType = typename expression::Result
<
    typename expression::AllContainers<Expression>::type
>::template elementType<Expression>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
