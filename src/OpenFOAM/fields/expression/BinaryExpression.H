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
    Foam::BinaryExpression

Description
    Binary expression class. Wraps a two sub-types which may be other
    expressions or expressionable containers or constant values. Applies a
    binary operation to the output of these sub-expressions.

    Note the specialisation of this class. This stores a sub-expression if it
    is provided as an r-value, thereby ensuring that it remains in scope as
    long as the outer expression does. And it takes a constant reference if a
    sub-expression is provided as a reference, which ensures that copies are
    not made unnecessarily.

\*---------------------------------------------------------------------------*/

#ifndef BinaryExpression_H
#define BinaryExpression_H

#include "expressionContainers.H"
#include "expressionOperate.H"
#include "expressionName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class BinaryExpression Declaration
\*---------------------------------------------------------------------------*/

template<class A, class B, class Op>
struct BinaryExpression
{
    A a;

    B b;

    BinaryExpression(A&& a, B&& b)
    :
        a(std::move(a)),
        b(std::move(b))
    {}
};

template<class A, class B, class Op>
class BinaryExpression<A, const B&, Op>
{
    const ConditionalIfBaseOf<refCount, B, tmp<B>, nil> tb_;

public:

    A a;

    const B& b;

    BinaryExpression(A&& a, const B& b)
    :
        a(std::move(a)),
        b(b)
    {}

    BinaryExpression(A&& a, const tmp<B>& tb)
    :
        tb_(tb, true),
        a(std::move(a)),
        b(tb_())
    {
        tb.clear();
    }
};

template<class A, class B, class Op>
class BinaryExpression<const A&, B, Op>
{
    const ConditionalIfBaseOf<refCount, A, tmp<A>, nil> ta_;

public:

    const A& a;

    B b;

    BinaryExpression(const A& a, B&& b)
    :
        a(a),
        b(std::move(b))
    {}

    BinaryExpression(const tmp<A>& ta, B&& b)
    :
        ta_(ta, true),
        a(ta_()),
        b(std::move(b))
    {
        ta.clear();
    }
};

template<class A, class B, class Op>
class BinaryExpression<const A&, const B&, Op>
{
    const ConditionalIfBaseOf<refCount, A, tmp<A>, nil> ta_;

    const ConditionalIfBaseOf<refCount, B, tmp<B>, nil> tb_;

public:

    const A& a;

    const B& b;

    BinaryExpression(const A& a, const B& b)
    :
        a(a),
        b(b)
    {}

    BinaryExpression(const tmp<A>& ta, const B& b)
    :
        ta_(ta, true),
        a(ta_()),
        b(b)
    {
        ta.clear();
    }

    BinaryExpression(const A& a, const tmp<B>& tb)
    :
        tb_(tb, true),
        a(a),
        b(tb_())
    {
        tb.clear();
    }

    BinaryExpression(const tmp<A>& ta, const tmp<B>& tb)
    :
        ta_(ta, true),
        tb_(tb, true),
        a(ta_()),
        b(tb_())
    {
        ta.clear();
        tb.clear();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace expression
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Helper type definitions and functions for forwarding the sub-expressions as
//  arguments. If the type is a binary expression, then forward the
//  sub-expression, otherwise forward the type itself.

template<class Type>
struct ArgA
{
    using type = Type;
};

template<class Type>
struct ArgB
{
    using type = Type;
};

template<class A, class B, class Op>
struct ArgA<BinaryExpression<A, B, Op>>
{
    using type = A;
};

template<class A, class B, class Op>
struct ArgB<BinaryExpression<A, B, Op>>
{
    using type = B;
};

template<class Type>
using ArgAType = typename ArgA<Type>::type;

template<class Type>
using ArgBType = typename ArgB<Type>::type;

template<class Type>
const Type& argA(const Type& t)
{
    return t;
}

template<class Type>
const Type& argB(const Type& t)
{
    return t;
}

template<class Type>
Type& argA(Type& t)
{
    return t;
}

template<class Type>
Type& argB(Type& t)
{
    return t;
}

template<class A, class B, class Op>
const A& argA(const BinaryExpression<A, B, Op>& e)
{
    return e.a;
}

template<class A, class B, class Op>
const B& argB(const BinaryExpression<A, B, Op>& e)
{
    return e.b;
}

template<class A, class B, class Op>
A& argA(BinaryExpression<A, B, Op>& e)
{
    return e.a;
}

template<class A, class B, class Op>
B& argB(BinaryExpression<A, B, Op>& e)
{
    return e.b;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialise Expressionable for the binary expression. A binary expression is
//  considered expressionable if either of it's sub-expressions are
//  expressionable.

template<class A, class B, class Op>
struct Expressionable<BinaryExpression<A, B, Op>>
:
    public std::conditional
    <
        (
            Expressionable<std::decay_t<A>>::value < 0
         || Expressionable<std::decay_t<B>>::value < 0
        ),
        notExpressionable,
        typename std::conditional
        <
            (
                Expressionable<std::decay_t<A>>::value
              > Expressionable<std::decay_t<B>>::value
            ),
            Expressionable<std::decay_t<A>>,
            Expressionable<std::decay_t<B>>
        >::type
    >::type
{};

template<class A, class B, class Op, class ... Args>
struct Expressionable<BinaryExpression<A, B, Op>, Arguments<Args ...>>
:
    public std::conditional
    <
        (
            Expressionable
            <
                std::decay_t<A>,
                Arguments<ArgAType<Args> ...>
            >::value
          > Expressionable
            <
                std::decay_t<B>,
                Arguments<ArgBType<Args> ...>
            >::value
        ),
        Expressionable
        <
            std::decay_t<A>,
            Arguments<ArgAType<Args> ...>
        >,
        Expressionable
        <
            std::decay_t<B>,
            Arguments<ArgBType<Args> ...>
        >
    >::type
{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialise AllContainers and TypeContainers for the binary expression so
//  that access is forwarded to the sub-expressions

template<class A, class B, class Op>
struct AllContainers
<
    BinaryExpression<A, B, Op>,
    EnableIfExpressionable<BinaryExpression<A, B, Op>>
>
{
    using allContainersA = AllContainers<std::decay_t<A>>;
    using allContainersB = AllContainers<std::decay_t<B>>;

    static const unsigned n = allContainersA::n + allContainersB::n;

    using type =
        typename std::conditional
        <
            (
                Expressionable<typename allContainersA::type>::value
              > Expressionable<typename allContainersB::type>::value
            ),
            typename allContainersA::type,
            typename allContainersB::type
        >::type;

    template<class Property, unsigned I>
    static std::enable_if_t
    <
        (I < allContainersA::n),
        typename Property::type
    > get
    (
        const BinaryExpression<A, B, Op>& e
    )
    {
        return allContainersA::template get<Property, I>(e.a);
    }

    template<class Property, unsigned I>
    static std::enable_if_t
    <
        (I >= allContainersA::n),
        typename Property::type
    > get
    (
        const BinaryExpression<A, B, Op>& e
    )
    {
        return allContainersB::template get<Property, I - allContainersA::n>
        (
            e.b
        );
    }

    template<class Property, unsigned I, unsigned N>
    static void get
    (
        const BinaryExpression<A, B, Op>& e,
        FixedList<typename Property::type, N>& result
    )
    {
        allContainersA::template get<Property, I>(e.a, result);
        allContainersB::template get<Property, allContainersA::n + I>
        (
            e.b,
            result
        );
    }
};

template<class A, class B, class Op, class Type>
struct TypeContainers
<
    BinaryExpression<A, B, Op>,
    Type,
    EnableIfNotBaseOf<Type, BinaryExpression<A, B, Op>>
>
{
    using typeContainersA = TypeContainers<std::decay_t<A>, Type>;
    using typeContainersB = TypeContainers<std::decay_t<B>, Type>;

    static const unsigned n = typeContainersA::n + typeContainersB::n;

    template<class Property, unsigned I>
    static std::enable_if_t
    <
        (I < typeContainersA::n),
        typename Property::type
    > get
    (
        const BinaryExpression<A, B, Op>& e
    )
    {
        return typeContainersA::template get<Property, I>(e.a);
    }

    template<class Property, unsigned I>
    static std::enable_if_t
    <
        (I >= typeContainersA::n),
        typename Property::type
    > get
    (
        const BinaryExpression<A, B, Op>& e
    )
    {
        return typeContainersB::template get<Property, I - typeContainersA::n>
        (
            e.b
        );
    }

    template<class Property, unsigned I, unsigned N>
    static void get
    (
        const BinaryExpression<A, B, Op>& e,
        FixedList<typename Property::type, N>& result
    )
    {
        typeContainersA::template get<Property, I>(e.a, result);
        typeContainersB::template get<Property, typeContainersA::n + I>
        (
            e.b,
            result
        );
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialise Operate for the binary expression so that operations are
//  forwarded to the sub-expressions.
//
//  Note different specialisations for expressionable and non-expressionable
//  binary expressions. The non-expressionable specialisation access returns a
//  value as calculated by the operation, whilst the expressionable
//  specialisation returns another expression containing the accessed
//  sub-types. In this way, the expressions can be nested, with each level
//  generating an appropriate type to which an element of an expressionable
//  container (container of containers, etc...) can be assigned.

template<class A, class B, class Op, class ... Args>
struct BinaryExpressionOperateBase
{
    using operateA =
        Operate<std::decay_t<A>, Arguments<ArgAType<Args> ...>>;
    using operateB =
        Operate<std::decay_t<B>, Arguments<ArgBType<Args> ...>>;

    struct compareA
    {
        static bool compare
        (
            const BinaryExpression<A, B, Op>& e,
            const Args& ... args
        )
        {
            return operateA::compare(e.a, argA(args) ...);
        }
    };

    struct compareB
    {
        static bool compare
        (
            const BinaryExpression<A, B, Op>& e,
            const Args& ... args
        )
        {
            return operateB::compare(e.b, argB(args) ...);
        }
    };

    static bool compare
    (
        const BinaryExpression<A, B, Op>& e,
        const Args& ... args
    )
    {
        return
            std::conditional
            <
                (Expressionable<std::decay_t<A>>::value > 0),
                compareA,
                compareB
            >::type::compare(e, args ...);
    }

    static void advance
    (
        const BinaryExpression<A, B, Op>& e,
        Args& ... args
    )
    {
        operateA::advance(e.a, argA(args) ...);
        operateB::advance(e.b, argB(args) ...);
    }
};

template<class A, class B, class Op, class ... Args>
struct Operate
<
    BinaryExpression<A, B, Op>,
    Arguments<Args ...>,
    EnableIfNotExpressionable<BinaryExpression<A, B, Op>, Args ...>
>
:
    public BinaryExpressionOperateBase<A, B, Op, Args ...>
{
    using operateA =
        typename BinaryExpressionOperateBase<A, B, Op, Args ...>::operateA;
    using operateB =
        typename BinaryExpressionOperateBase<A, B, Op, Args ...>::operateB;

    using type =
        decltype
        (
            Op::value
            (
                std::declval<typename operateA::type>(),
                std::declval<typename operateB::type>()
            )
        );

    static type access
    (
        const BinaryExpression<A, B, Op>& e,
        const Args& ... args
    )
    {
        return
            Op::value
            (
                operateA::access(e.a, argA(args) ...),
                operateB::access(e.b, argB(args) ...)
            );
    }
};

template<class A, class B, class Op, class ... Args>
struct Operate
<
    BinaryExpression<A, B, Op>,
    Arguments<Args ...>,
    EnableIfExpressionable<BinaryExpression<A, B, Op>, Args ...>
>
:
    public BinaryExpressionOperateBase<A, B, Op, Args ...>
{
    using operateA =
        typename BinaryExpressionOperateBase<A, B, Op, Args ...>::operateA;
    using operateB =
        typename BinaryExpressionOperateBase<A, B, Op, Args ...>::operateB;

    using type =
        BinaryExpression
        <
            typename operateA::type,
            typename operateB::type,
            Op
        >;

    static type access
    (
        const BinaryExpression<A, B, Op>& e,
        const Args& ... args
    )
    {
        return
            type
            (
                operateA::access(e.a, argA(args) ...),
                operateB::access(e.b, argB(args) ...)
            );
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialise Operation for the binary expression. Take the base operator type
//  and precedence from the operation parameter.

template<class A, class B, class Op>
struct Operation<BinaryExpression<A, B, Op>>
{
    using type = Op;

    static const signed char precedence = Op::precedence;
};

//- Specialise Name for the binary expression. Wraps the names of the
//  sub-expressions in the name of the operation, adding brackets if the
//  relative precedence of the operations requires it.

template<class A, class B, class Op>
struct Name
<
    BinaryExpression<A, B, Op>,
    EnableIfExpressionable<BinaryExpression<A, B, Op>>
>
{
    using nameA = Name<std::decay_t<A>>;
    using nameB = Name<std::decay_t<B>>;

    static word value(const BinaryExpression<A, B, Op>& e)
    {
        return
            Op::name
            (
                Op::precedence != Operation<void>::precedence
             && Operation<A>::precedence != Operation<void>::precedence
             && Op::precedence > Operation<A>::precedence
              ? '(' + nameA::value(e.a) + ')'
              : nameA::value(e.a),
                Op::precedence != Operation<void>::precedence
             && Operation<B>::precedence != Operation<void>::precedence
             && (
                    Op::precedence > Operation<B>::precedence
                 || (
                        Op::precedence == Operation<B>::precedence
                     && !Op::associative
                    )
                )
              ? '(' + nameB::value(e.b) + ')'
              : nameB::value(e.b)
            );
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
