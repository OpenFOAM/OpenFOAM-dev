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
    Unary and binary scalar functions and operators used by the expression
    templating system

\*---------------------------------------------------------------------------*/

#ifndef scalarExpressions_H
#define scalarExpressions_H

#include "scalar.H"
#include "expressions.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class powOp Declaration
\*---------------------------------------------------------------------------*/

struct powOp
{
    static inline const char* prefix = "pow(";
    static inline const char infix = ',';
    static inline const char postfix = ')';

    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static scalar value(const scalar a, const float b)
    {
        return pow(a, b);
    }

    static scalar value(const scalar a, const double b)
    {
        return pow(a, b);
    }

    static scalar value(const scalar a, const long double b)
    {
        return pow(a, b);
    }

    static scalar value(const scalar a, const int32_t b)
    {
        return integerPow(a, b);
    }

    static scalar value(const scalar a, const int64_t b)
    {
        return integerPow(a, b);
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return pow(a, b);
    }
};

/*---------------------------------------------------------------------------*\
                        Class integerPowOp Declaration
\*---------------------------------------------------------------------------*/

struct integerPowOp
{
    static inline const char* prefix = "integerPow(";
    static inline const char infix = ',';
    static inline const char postfix = ')';

    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static scalar value(const scalar a, const label b)
    {
        return integerPow(a, b);
    }
};

/*---------------------------------------------------------------------------*\
                        Class integerRootOp Declaration
\*---------------------------------------------------------------------------*/

struct integerRootOp
{
    static inline const char* prefix = "integerRoot(";
    static inline const char infix = ',';
    static inline const char postfix = ')';

    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static scalar value(const scalar a, const label b)
    {
        return integerRoot(a, b);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Unary scalar functions

UnaryOpFunction(pow3)
UnaryOpFunction(pow4)
UnaryOpFunction(pow5)
UnaryOpFunction(pow6)
UnaryOpFunction(pow025)
UnaryOpFunction(sqrt)
UnaryOpFunction(sign)
UnaryOpFunction(pos)
UnaryOpFunction(pos0)
UnaryOpFunction(neg)
UnaryOpFunction(neg0)
UnaryOpFunction(posPart)
UnaryOpFunction(negPart)
UnaryOpFunction(exp)
UnaryOpFunction(log)
UnaryOpFunction(log10)
UnaryOpFunction(sin)
UnaryOpFunction(cos)
UnaryOpFunction(tan)
UnaryOpFunction(asin)
UnaryOpFunction(acos)
UnaryOpFunction(atan)

// Expensive unary scalar functions
UnaryOpFunction(cbrt)
UnaryOpFunction(sinh)
UnaryOpFunction(cosh)
UnaryOpFunction(tanh)
UnaryOpFunction(asinh)
UnaryOpFunction(acosh)
UnaryOpFunction(atanh)
UnaryOpFunction(erf)
UnaryOpFunction(erfc)
UnaryOpFunction(lgamma)
UnaryOpFunction(j0)
UnaryOpFunction(j1)
UnaryOpFunction(y0)
UnaryOpFunction(y1)

BinaryFunction(stabilise);
BinaryFunction(integerPow);
BinaryFunction(integerRoot);

// Eager evaluate pow to avoid interference with loop optimisation
// Benefit demonstrated on x86 processors with sse2 and avx512
BinaryEagerFunction(pow);

/*---------------------------------------------------------------------------*\
                            Class atan2Op Declaration
\*---------------------------------------------------------------------------*/

struct atan2Op
{
    static inline const char* prefix = "atan2(";
    static inline const char infix = ',';
    static inline const char postfix = ')';

    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static scalar value(const scalar a, const scalar b)
    {
        return atan2(a, b);
    }
};

// Atan2 function
BinaryFunction(atan2);

/*---------------------------------------------------------------------------*\
                            Class jnOp Declaration
\*---------------------------------------------------------------------------*/

struct jnOp
{
    static inline const char* prefix = "jn(";
    static inline const char infix = ',';
    static inline const char postfix = ')';

    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static scalar value(const int n, const scalar s)
    {
        return jn(n, s);
    }
};

// Bessel function jn
BinaryFunction(jn);

/*---------------------------------------------------------------------------*\
                            Class ynOp Declaration
\*---------------------------------------------------------------------------*/

struct ynOp
{
    static inline const char* prefix = "yn(";
    static inline const char infix = ',';
    static inline const char postfix = ')';

    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static scalar value(const int n, const scalar s)
    {
        return yn(n, s);
    }
};

// Bessel function yn
BinaryFunction(yn);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
