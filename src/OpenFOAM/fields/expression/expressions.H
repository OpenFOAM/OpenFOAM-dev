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

Description
    Unary and binary functions and operators used by the expression
    templating system

\*---------------------------------------------------------------------------*/

#ifndef expressions_H
#define expressions_H

#include "expressionFunctor.H"
#include "expressionReduce.H"
#include "PstreamReduceOps.H"
#include "expressionEval.H"

#include "UnaryExpression.H"
#include "BinaryExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define a unary function or operator given the global function or operator
// name and the operator class that it uses

#define UnaryFunction__(Func, Op, XPar, XArg, XVar)                            \
    template                                                                   \
    <                                                                          \
        class X,                                                               \
        class = EnableIfExpressionable<UnaryExpression<X, Op>>                 \
    >                                                                          \
    UnaryExpression<XPar, Op> Func(XArg)                                       \
    {                                                                          \
        return UnaryExpression<XPar, Op>(XVar);                                \
    }                                                                          \

#define UnaryFunction_(Func, Op)                                               \
    UnaryFunction__(Func, Op, const X&, const X& x, x)                         \
    UnaryFunction__(Func, Op, const X&, X& x, x)                               \
    UnaryFunction__(Func, Op, X, X&& x, std::move(x))                          \
    UnaryFunction__(Func, Op, const X&, const tmp<X>& x, x)                    \
    UnaryFunction__(Func, Op, const X&, tmp<X>& x, x)                          \
    UnaryFunction__(Func, Op, const X&, tmp<X>&& x, x)

#define UnaryOperator(OpSym, Func)                                             \
    UnaryFunction_(operator OpSym, Func##Op)

#define UnaryFunction(Func)                                                    \
    UnaryFunction_(Func, Func##Op)

#define UnaryOpFunction(RetType, Type, Func)                                   \
                                                                               \
struct Func##Op                                                                \
{                                                                              \
    static const signed char precedence =                                      \
        expression::Operation<void>::precedence;                               \
                                                                               \
    static word name(const word& x)                                            \
    {                                                                          \
        return #Func"(" + x + ')';                                             \
    }                                                                          \
                                                                               \
    static RetType value(const Type& x)                                        \
    {                                                                          \
        return Func(x);                                                        \
    }                                                                          \
                                                                               \
    static dimensionSet value(const dimensionSet& x)                           \
    {                                                                          \
        return Func(x);                                                        \
    }                                                                          \
};                                                                             \
                                                                               \
UnaryFunction(Func)

// Promotion operator (i.e., unary plus). Has no effect. Just facilitates
// expressions with plus signs added for clarity.
UnaryOperator(+, promote)

// Negation operator (i.e., unary minus). Changes the sign of the operand.
UnaryOperator(-, negate)

UnaryFunction(sqr)
UnaryFunction(mag)
UnaryFunction(magSqr)
UnaryFunction(cmptMax)
UnaryFunction(cmptMin)
UnaryFunction(cmptAv)
UnaryFunction(cmptMag)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define a unary reduction function given the global function name and the
// (binary) operator class that it uses

#define UnaryReduction_(Func, gFunc, Op)                                       \
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
    >::template elementType<Expression> Func(const tmp<Expression>& te)        \
    {                                                                          \
        auto x = Func(te());                                                   \
        te.clear();                                                            \
        return x;                                                              \
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
    }                                                                          \
                                                                               \
    template<class Expression, class = EnableIfExpressionable<Expression>>     \
    typename expression::Result                                                \
    <                                                                          \
        typename expression::AllContainers<Expression>::type                   \
    >::template elementType<Expression> gFunc(const tmp<Expression>& te)       \
    {                                                                          \
        auto x = gFunc(te());                                                  \
        te.clear();                                                            \
        return x;                                                              \
    }

#define UnaryReduction(Func, gFunc)                                            \
    UnaryReduction_(Func, gFunc, Func##Op)

UnaryReduction(sum, gSum);
UnaryReduction(max, gMax);
UnaryReduction(min, gMin);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define a binary function or operator given the global function or operator
// name and the operator class that it uses

#define BinaryFunction___(Func, Op, APar, AArg, AVar, BPar, BArg, BVar)        \
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

#define BinaryFunction__(Func, Op, APar, AArg, AVar)                           \
    BinaryFunction___(Func, Op, APar, AArg, AVar, const B&, const B& b, b)     \
    BinaryFunction___(Func, Op, APar, AArg, AVar, const B&, B& b, b)           \
    BinaryFunction___(Func, Op, APar, AArg, AVar, B, B&& b, std::move(b))      \
    BinaryFunction___(Func, Op, APar, AArg, AVar, const B&, const tmp<B>& b, b)\
    BinaryFunction___(Func, Op, APar, AArg, AVar, const B&, tmp<B>& b, b)      \
    BinaryFunction___(Func, Op, APar, AArg, AVar, const B&, tmp<B>&& b, b)

#define BinaryFunction_(Func, Op)                                              \
    BinaryFunction__(Func, Op, const A&, const A& a, a)                        \
    BinaryFunction__(Func, Op, const A&, A& a, a)                              \
    BinaryFunction__(Func, Op, A, A&& a, std::move(a))                         \
    BinaryFunction__(Func, Op, const A&, const tmp<A>& a, a)                   \
    BinaryFunction__(Func, Op, const A&, tmp<A>& a, a)                         \
    BinaryFunction__(Func, Op, const A&, tmp<A>&& a, a)

#define BinaryOperator(OpSym, Func)                                            \
    BinaryFunction_(operator OpSym, Func##Op)

#define BinaryFunction(Func)                                                   \
    BinaryFunction_(Func, Func##Op)

// Addition operator and sum function
BinaryOperator(+, add);

// Subtraction operator
BinaryOperator(-, subtract);

// Multiplication and/or outer-product operator
BinaryOperator(*, multiply);

// Division operator
BinaryOperator(/, divide);

// Inner product (dot) operator
BinaryOperator(&, inner);

// Cross product operator
BinaryOperator(^, cross);

// Maximum functions
BinaryFunction(max);

// Minimum function
BinaryFunction(min);

BinaryFunction(cmptMultiply);
BinaryFunction(cmptDivide);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
