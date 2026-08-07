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
    Foam::Op*

Description
    Operation classes, functions and operators used by the expression
    templating system

\*---------------------------------------------------------------------------*/

#ifndef operations_H
#define operations_H

#include "PstreamReduceOps.H"

#include "expressionFunctor.H"
#include "expressionReduce.H"

#include "UnaryExpression.H"
#include "BinaryExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define a unary function or operator given the global function or operator
//  name and the operator class that it uses

#define UnaryFunction_(Func, Op, XPar, XArg, XVar)                             \
    template                                                                   \
    <                                                                          \
        class X,                                                               \
        class = EnableIfExpressionable<UnaryExpression<X, Op>>                 \
    >                                                                          \
    UnaryExpression<XPar, Op> Func(XArg)                                       \
    {                                                                          \
        return UnaryExpression<XPar, Op>(XVar);                                \
    }                                                                          \

#define UnaryFunction(Func, Op)                                                \
    UnaryFunction_(Func, Op, const X&, const X& x, x)                          \
    UnaryFunction_(Func, Op, const X&, X& x, x)                                \
    UnaryFunction_(Func, Op, X, X&& x, std::move(x))                           \
    UnaryFunction_(Func, Op, const X&, const tmp<X>& x, x)                     \
    UnaryFunction_(Func, Op, const X&, tmp<X>& x, x)                           \
    UnaryFunction_(Func, Op, const X&, tmp<X>&& x, x)

/*---------------------------------------------------------------------------*\
                          Class OpPromote Declaration
\*---------------------------------------------------------------------------*/

struct OpPromote
{
    static const char precedence = 2;

    static word name(const word& x)
    {
        return x;
    }

    template<class X>
    static X value(const X& x)
    {
        return x;
    }
};

//- Promotion operator (i.e., unary plus). Has no effect. Just facilitates
//  expressions with plus signs added for clarity.
UnaryFunction(operator+, OpPromote)

/*---------------------------------------------------------------------------*\
                          Class OpNegate Declaration
\*---------------------------------------------------------------------------*/

struct OpNegate
{
    static const char precedence = 2;

    static word name(const word& x)
    {
        return '-' + x;
    }

    template<class X>
    static auto value(const X& x)
    {
        return - x;
    }
};

//- Negation operator (i.e., unary minus). Changes the sign of the operand.
UnaryFunction(operator-, OpNegate)

/*---------------------------------------------------------------------------*\
                            Class OpSqr Declaration
\*---------------------------------------------------------------------------*/

struct OpSqr
{
    static const char precedence =
        expression::Operation<void>::precedence;

    static word name(const word& x)
    {
        return "sqr(" + x + ')';
    }

    template<class X>
    static auto value(const X& x)
    {
        return x*x;
    }
};

//- Square function
UnaryFunction(sqr, OpSqr)

/*---------------------------------------------------------------------------*\
                            Class OpSqrt Declaration
\*---------------------------------------------------------------------------*/

struct OpSqrt
{
    static const char precedence =
        expression::Operation<void>::precedence;

    static word name(const word& x)
    {
        return "sqrt(" + x + ')';
    }

    static scalar value(const scalar x)
    {
        return sqrt(x);
    }
};

//- Square root function
UnaryFunction(sqrt, OpSqrt)

/*---------------------------------------------------------------------------*\
                            Class OpMag Declaration
\*---------------------------------------------------------------------------*/

struct OpMag
{
    static const char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& x)
    {
        return "mag(" + x + ')';
    }

    template<class X>
    static auto value(const X& x)
    {
        return mag(x);
    }
};

//- Magnitude function
UnaryFunction(mag, OpMag)

/*---------------------------------------------------------------------------*\
                          Class OpMagSqr Declaration
\*---------------------------------------------------------------------------*/

struct OpMagSqr
{
    static const char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& x)
    {
        return "magSqr(" + x + ')';
    }

    template<class X>
    static auto value(const X& x)
    {
        return magSqr(x);
    }
};

//- Magnitude-squared function
UnaryFunction(magSqr, OpMagSqr)

/*---------------------------------------------------------------------------*\
                            Class OpLog Declaration
\*---------------------------------------------------------------------------*/

struct OpLog
{
    static const char precedence =
        expression::Operation<void>::precedence;

    static word name(const word& x)
    {
        return "log(" + x + ')';
    }

    static scalar value(const scalar x)
    {
        return log(x);
    }
};

//- Natural logarithm function
UnaryFunction(log, OpLog)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef UnaryFunction_
#undef UnaryFunction

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define a binary function or operator given the global function or operator
//  name and the operator class that it uses

#define BinaryFunction__(Func, Op, APar, AArg, AVar, BPar, BArg, BVar)         \
    template                                                                   \
    <                                                                          \
        class A,                                                               \
        class B,                                                               \
        class = EnableIfExpressionable<BinaryExpression<A, B, Op>>             \
    >                                                                          \
    BinaryExpression<APar, BPar, Op> Func(AArg, BArg)                          \
    {                                                                          \
        return BinaryExpression<APar, BPar, Op>(AVar, BVar);                   \
    }                                                                          \

#define BinaryFunction_(Func, Op, APar, AArg, AVar)                            \
    BinaryFunction__(Func, Op, APar, AArg, AVar, const B&, const B& b, b)      \
    BinaryFunction__(Func, Op, APar, AArg, AVar, const B&, B& b, b)            \
    BinaryFunction__(Func, Op, APar, AArg, AVar, B, B&& b, std::move(b))       \
    BinaryFunction__(Func, Op, APar, AArg, AVar, const B&, const tmp<B>& b, b) \
    BinaryFunction__(Func, Op, APar, AArg, AVar, const B&, tmp<B>& b, b)       \
    BinaryFunction__(Func, Op, APar, AArg, AVar, const B&, tmp<B>&& b, b)

#define BinaryFunction(Func, Op)                                               \
    BinaryFunction_(Func, Op, const A&, const A& a, a)                         \
    BinaryFunction_(Func, Op, const A&, A& a, a)                               \
    BinaryFunction_(Func, Op, A, A&& a, std::move(a))                          \
    BinaryFunction_(Func, Op, const A&, const tmp<A>& a, a)                    \
    BinaryFunction_(Func, Op, const A&, tmp<A>& a, a)                          \
    BinaryFunction_(Func, Op, const A&, tmp<A>&& a, a)

//- Define a unary reduction function given the global function name and the
//  (binary) operator class that it uses

#define UnaryReduction(Func, gFunc, Op)                                        \
                                                                               \
    template<class Expression, class = EnableIfExpressionable<Expression>>     \
    typename expression::Result                                                \
    <                                                                          \
        typename expression::AllContainers<Expression>::type                   \
    >::template elementType<Expression> Func(const Expression& e)              \
    {                                                                          \
        return                                                                 \
            expression::Reduce                                                 \
            <                                                                  \
                typename expression::AllContainers<Expression>::type           \
            >::template value<Op>(e);                                          \
    }                                                                          \
                                                                               \
    template<class Expression, class = EnableIfExpressionable<Expression>>     \
    typename expression::Result                                                \
    <                                                                          \
        typename expression::AllContainers<Expression>::type                   \
    >::template elementType<Expression> gFunc                                  \
    (                                                                          \
        const Expression& e,                                                   \
        const label comm = UPstream::worldComm                                 \
    )                                                                          \
    {                                                                          \
        /* local reduction */                                                  \
        auto x =                                                               \
            expression::Reduce                                                 \
            <                                                                  \
                typename expression::AllContainers<Expression>::type           \
            >::template value<Op>(e);                                          \
                                                                               \
        /* global reduction */                                                 \
        reduce(x, expression::Functor<Op>(), Pstream::msgType(), comm);        \
                                                                               \
        return x;                                                              \
    }

/*---------------------------------------------------------------------------*\
                            Class OpAdd Declaration
\*---------------------------------------------------------------------------*/

struct OpAdd
{
    static const char precedence = 0;

    static const bool associative = true;

    static word name(const word& a, const word& b)
    {
        return a + '+' + b;
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return a + b;
    }

    static word name(const word& x)
    {
        return "sum(" + x + ')';
    }

    template<class X>
    static const X& identity()
    {
        return pTraits<X>::zero;
    }
};

//- Addition operator and sum function
BinaryFunction(operator+, OpAdd);
UnaryReduction(sum, gSum, OpAdd);

/*---------------------------------------------------------------------------*\
                          Class OpSubtract Declaration
\*---------------------------------------------------------------------------*/

struct OpSubtract
{
    static const char precedence = 0;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return a + '-' + b;
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return a - b;
    }
};

//- Subtraction operator
BinaryFunction(operator-, OpSubtract);

/*---------------------------------------------------------------------------*\
                          Class OpMultiply Declaration
\*---------------------------------------------------------------------------*/

struct OpMultiply
{
    static const char precedence = 1;

    static const bool associative = true;

    static word name(const word& a, const word& b)
    {
        return a + '*' + b;
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return a*b;
    }
};

//- Multiplication and/or outer-product operator
BinaryFunction(operator*, OpMultiply);

/*---------------------------------------------------------------------------*\
                            Class OpDivide Declaration
\*---------------------------------------------------------------------------*/

struct OpDivide
{
    static const char precedence = 1;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return a + '|' + b;
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return a/b;
    }
};

//- Division operator
BinaryFunction(operator/, OpDivide);

/*---------------------------------------------------------------------------*\
                            Class OpInner Declaration
\*---------------------------------------------------------------------------*/

struct OpInner
{
    static const char precedence = 1;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return a + '&' + b;
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return a & b;
    }
};

//- Inner product (dot) operator
BinaryFunction(operator&, OpInner);

/*---------------------------------------------------------------------------*\
                            Class OpCross Declaration
\*---------------------------------------------------------------------------*/

struct OpCross
{
    static const char precedence = 1;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return a + '^' + b;
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return a ^ b;
    }
};

//- Cross product operator
BinaryFunction(operator^, OpCross);

/*---------------------------------------------------------------------------*\
                            Class OpPower Declaration
\*---------------------------------------------------------------------------*/

struct OpPower
{
    static const char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return "pow(" + a + ',' + b + ')';
    }

    static scalar value(const scalar a, const scalar b)
    {
        return pow(a, b);
    }

    static scalar value(const scalar a, const label b)
    {
        if (b == 0) return 1;

        scalar aa = b > 0 ? a : 1/a, result = 1;

        for (label i = b; i != 0; i /= 2)
        {
            if (i % 2) result *= aa;

            aa *= aa;
        }

        return result;
    }
};

//- Power function
BinaryFunction(pow, OpPower);

/*---------------------------------------------------------------------------*\
                             Class OpMax Declaration
\*---------------------------------------------------------------------------*/

struct OpMax
{
    static const char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return "max(" + a + ',' + b + ')';
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return max(a, b);
    }

    static word name(const word& x)
    {
        return "max(" + x + ')';
    }

    template<class X>
    static const X& identity()
    {
        return pTraits<X>::min;
    }
};

//- Maximum functions
BinaryFunction(max, OpMax);
UnaryReduction(max, gMax, OpMax);

/*---------------------------------------------------------------------------*\
                             Class OpMin Declaration
\*---------------------------------------------------------------------------*/

struct OpMin
{
    static const char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return "min(" + a + ',' + b + ')';
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return min(a, b);
    }

    static word name(const word& x)
    {
        return "min(" + x + ')';
    }

    template<class X>
    static const X& identity()
    {
        return pTraits<X>::max;
    }
};

//- Minimum function
BinaryFunction(min, OpMin);
UnaryReduction(min, gMin, OpMin);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef BinaryFunction__
#undef BinaryFunction_
#undef BinaryFunction
#undef UnaryReduction

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
