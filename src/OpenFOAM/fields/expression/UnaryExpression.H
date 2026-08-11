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
    Foam::UnaryExpression

Description
    Unary expression class. Wraps a single sub-type which may be another
    expression or an expressionable container or a constant value. Applies a
    unary operation to the output of this sub-expression.

    Note the specialisation of this class. This stores the sub-expression if it
    is provided as an r-value, thereby ensuring that it remains in scope as
    long as the outer expression does. And it takes a constant reference if the
    sub-expression is provided as a reference, which ensures that copies are
    not made unnecessarily.

\*---------------------------------------------------------------------------*/

#ifndef UnaryExpression_H
#define UnaryExpression_H

#include "expressionContainers.H"
#include "expressionOperate.H"
#include "expressionName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class UnaryExpression Declaration
\*---------------------------------------------------------------------------*/

template<class X, class Op>
struct UnaryExpression
{
    X x;

    UnaryExpression(X&& x)
    :
        x(std::move(x))
    {}
};

template<class X, class Op>
class UnaryExpression<const X&, Op>
{
    const ConditionalIfBaseOf<refCount, X, tmp<X>, nil> tx_;

public:

    const X& x;

    UnaryExpression(const X& x)
    :
        x(x)
    {}

    UnaryExpression(const tmp<X>& tx)
    :
        tx_(tx, true),
        x(tx_())
    {
        tx.clear();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace expression
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Helper type definitions and functions for forwarding the sub-expression as
//  an argument. If the type is a unary expression, then forward the
//  sub-expression, otherwise forward the type itself.

template<class Type>
struct ArgX
{
    using type = Type;
};

template<class X, class Op>
struct ArgX<UnaryExpression<X, Op>>
{
    using type = X;
};

template<class Type>
using ArgXType = typename ArgX<Type>::type;

template<class Type>
const Type& argX(const Type& t)
{
    return t;
}

template<class Type>
Type& argX(Type& t)
{
    return t;
}

template<class X, class Op>
const X& argX(const UnaryExpression<X, Op>& e)
{
    return e.x;
}

template<class X, class Op>
X& argX(UnaryExpression<X, Op>& e)
{
    return e.x;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialise Expressionable for the unary expression. A unary expression is
//  considered expressionable if it's sub-expression is expressionable.

template<class X, class Op>
struct Expressionable<UnaryExpression<X, Op>>
:
    public Expressionable<std::decay_t<X>>
{};

template<class X, class Op, class ... Args>
struct Expressionable<UnaryExpression<X, Op>, Arguments<Args ...>>
:
    public Expressionable<std::decay_t<X>, Arguments<ArgXType<Args> ...>>
{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialise AllContainers and TypeContainers for the unary expression so
//  that access is forwarded to the sub-expression

template<class X, class Op>
struct AllContainers
<
    UnaryExpression<X, Op>,
    EnableIfExpressionable<UnaryExpression<X, Op>>
>
{
    using allContainersX = AllContainers<std::decay_t<X>>;

    static const unsigned n = allContainersX::n;

    using type = typename allContainersX::type;

    template<class Property, unsigned I>
    static typename Property::type get
    (
        const UnaryExpression<X, Op>& e
    )
    {
        return allContainersX::template get<Property, I>(e.x);
    }

    template<class Property, unsigned I, unsigned N>
    static void get
    (
        const UnaryExpression<X, Op>& e,
        FixedList<typename Property::type, N>& result
    )
    {
        allContainersX::template get<Property, I>(e.x, result);
    }
};

template<class X, class Op, class Type>
struct TypeContainers
<
    UnaryExpression<X, Op>,
    Type,
    EnableIfNotBaseOf<Type, UnaryExpression<X, Op>>
>
{
    using typeContainersX = TypeContainers<std::decay_t<X>, Type>;

    static const unsigned n = typeContainersX::n;

    template<class Property, unsigned I>
    static typename Property::type get
    (
        const UnaryExpression<X, Op>& e
    )
    {
        return typeContainersX::template get<Property, I>(e.x);
    }

    template<class Property, unsigned I, unsigned N>
    static void get
    (
        const UnaryExpression<X, Op>& e,
        FixedList<typename Property::type, N>& result
    )
    {
        typeContainersX::template get<Property, I>(e.x, result);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialise Operate for the unary expression so that operations are
//  forwarded to the sub-expression.
//
//  Note different specialisations for expressionable and non-expressionable
//  unary expressions. The non-expressionable specialisation access returns a
//  value as calculated by the operation, whilst the expressionable
//  specialisation returns another expression containing the accessed
//  sub-types. In this way, the expressions can be nested, with each level
//  generating an appropriate type to which an element of an expressionable
//  container (container of containers, etc...) can be assigned.

template<class X, class Op, class ... Args>
struct UnaryExpressionOperateBase
{
    using operateX =
        Operate<std::decay_t<X>, Arguments<ArgXType<Args> ...>>;

    static bool compare
    (
        const UnaryExpression<X, Op>& e,
        const Args& ... args
    )
    {
        return operateX::compare(e.x, argX(args) ...);
    }

    static void advance
    (
        const UnaryExpression<X, Op>& e,
        Args& ... args
    )
    {
        operateX::advance(e.x, argX(args) ...);
    }
};

template<class X, class Op, class ... Args>
struct Operate
<
    UnaryExpression<X, Op>,
    Arguments<Args ...>,
    EnableIfNotExpressionable<UnaryExpression<X, Op>, Args ...>
>
:
    public UnaryExpressionOperateBase<X, Op, Args ...>
{
    using operateX =
        typename UnaryExpressionOperateBase<X, Op, Args ...>::operateX;

    using type =
        decltype(Op::value(std::declval<typename operateX::type>()));

    static type access
    (
        const UnaryExpression<X, Op>& e,
        const Args& ... args
    )
    {
        return Op::value(operateX::access(e.x, argX(args) ...));
    }
};

template<class X, class Op, class ... Args>
struct Operate
<
    UnaryExpression<X, Op>,
    Arguments<Args ...>,
    EnableIfExpressionable<UnaryExpression<X, Op>, Args ...>
>
:
    public UnaryExpressionOperateBase<X, Op, Args ...>
{
    using operateX =
        typename UnaryExpressionOperateBase<X, Op, Args ...>::operateX;

    using type = UnaryExpression<typename operateX::type, Op>;

    static type access
    (
        const UnaryExpression<X, Op>& e,
        const Args& ... args
    )
    {
        return type(operateX::access(e.x, argX(args) ...));
    }
};

//- Specialise Operation for the unary expression. Take the base operator type
//  and precedence from the operation parameter.

template<class X, class Op>
struct Operation<UnaryExpression<X, Op>>
{
    using type = Op;

    static const signed char precedence = Op::precedence;
};

//- Specialise Name for the unary expression. Wraps the name of the
//  sub-expression in the name of the operation, adding brackets if the
//  relative precedence of the operations requires it.

template<class X, class Op>
struct Name
<
    UnaryExpression<X, Op>,
    EnableIfExpressionable<UnaryExpression<X, Op>>
>
{
    using nameX = Name<std::decay_t<X>>;

    static word value(const UnaryExpression<X, Op>& e)
    {
        return
            Op::name
            (
                Op::precedence != Operation<void>::precedence
             && Operation<X>::precedence != Operation<void>::precedence
             && Op::precedence > Operation<X>::precedence
              ? '(' + nameX::value(e.x) + ')'
              : nameX::value(e.x)
            );
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
