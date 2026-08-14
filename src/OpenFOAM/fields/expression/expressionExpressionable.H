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

Class
    Foam::expression::Expressionable

Description
    Trait which defines whether a type operates in the manner of an expression.
    It should be specialised for expressionable containers. It also optionally
    takes as parameters the list of arguments relevant to the current
    operation. This is needed to determine whether an intermediate expression
    (e.g., an iterator) produces an expressionable result or not.

\*---------------------------------------------------------------------------*/

#ifndef expressionExpressionable_H
#define expressionExpressionable_H

#include <algorithm>
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

/*---------------------------------------------------------------------------*\
                          Class expression Declaration
\*---------------------------------------------------------------------------*/

//- Base class for all unary and binary expressions
//  Used by EnableIfExpression
struct expression
{};


/*---------------------------------------------------------------------------*\
                          Class Arguments Declaration
\*---------------------------------------------------------------------------*/

template<class ... Args>
struct Arguments
{};


/*---------------------------------------------------------------------------*\
                        Class Expressionable Declaration
\*---------------------------------------------------------------------------*/

template<class Expression, class = Arguments<>>
struct Expressionable;

template<class Expression>
struct Expressionable<Expression>
:
    public std::integral_constant<signed char, 0>
{};

template<class Expression, class ... Args>
struct Expressionable<Expression, Arguments<Args ...>>
:
    public Expressionable<Expression>
{};


//- Not expressionable
using notExpressionable = std::integral_constant<signed char, -1>;

//- Not expressionable if the Expression is not expressionable
template<class Expression>
using NotExpressionable =
    std::integral_constant
    <
        signed char,
        std::min(static_cast<signed char>(0), Expressionable<Expression>::value)
    >;


/*---------------------------------------------------------------------------*\
                     Class ExpressionableIterable Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class = Arguments<>>
struct ExpressionableIterable;

//- Iterable containers are one more level expressionable than their type
template<class Type>
struct ExpressionableIterable<Type, Arguments<>>
:
    public std::integral_constant
    <
        signed char,
        std::max
        (
            static_cast<signed char>(0),
            Expressionable<typename Type::value_type>::value
        )
      + static_cast<signed char>(1)
    >
{};

//- Access by index to an iterable container is expressionable if the type is
template<class Type>
struct ExpressionableIterable<Type, Arguments<label>>
:
    public Expressionable<typename Type::value_type>
{};

//- Access by iterator to an iterable container is expressionable if the type is
template<class Type>
struct ExpressionableIterable
<
    Type,
    Arguments<typename Type::const_iterator>
>
:
    public Expressionable<typename Type::value_type>
{};

//- Other arguments are as expressionable as the container itself
template<class Type, class ... Args>
struct ExpressionableIterable<Type, Arguments<Args ...>>
:
    public Expressionable<Type>
{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Alias to enable instantiation if a type is a unary or binary expression
template<class Expression>
using EnableIfExpression =
    std::enable_if_t
    <
        std::is_base_of<expression::expression, Expression>::value,
        void
    >;

//- Alias to enable instantiation if a type (and arguments) is an expression
template<class Expression, class ... Args>
using EnableIfExpressionable =
    std::enable_if_t
    <
        (
            expression::Expressionable
            <
                Expression,
                expression::Arguments<Args ...>
            >::value > 0
        )
    >;

//- Alias to enable instantiation if a type (and arguments) is not an expression
template<class Expression, class ... Args>
using EnableIfNotExpressionable =
    std::enable_if_t
    <
        (
            expression::Expressionable
            <
                Expression,
                expression::Arguments<Args ...>
            >::value <= 0
        )
    >;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
