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
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return "pow(" + a + ',' + b + ')';
    }

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

    template<class Exp>
    static dimensionSet value(const dimensionSet& a, const Exp& b)
    {
        return pow(a, b);
    }
};

/*---------------------------------------------------------------------------*\
                        Class integerPowOp Declaration
\*---------------------------------------------------------------------------*/

struct integerPowOp
{
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return "integerPow(" + a + ',' + b + ')';
    }

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
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return "integerRoot(" + a + ',' + b + ')';
    }

    static scalar value(const scalar a, const label b)
    {
        return integerRoot(a, b);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Unary scalar functions

#define UnaryScalarOpFunction(func)                                            \
    UnaryOpFunction(scalar, scalar, func)

UnaryScalarOpFunction(pow3)
UnaryScalarOpFunction(pow4)
UnaryScalarOpFunction(pow5)
UnaryScalarOpFunction(pow6)
UnaryScalarOpFunction(pow025)
UnaryScalarOpFunction(sqrt)
UnaryScalarOpFunction(cbrt)
UnaryScalarOpFunction(sign)
UnaryScalarOpFunction(pos)
UnaryScalarOpFunction(pos0)
UnaryScalarOpFunction(neg)
UnaryScalarOpFunction(neg0)
UnaryScalarOpFunction(posPart)
UnaryScalarOpFunction(negPart)
UnaryScalarOpFunction(exp)
UnaryScalarOpFunction(log)
UnaryScalarOpFunction(log10)
UnaryScalarOpFunction(sin)
UnaryScalarOpFunction(cos)
UnaryScalarOpFunction(tan)
UnaryScalarOpFunction(asin)
UnaryScalarOpFunction(acos)
UnaryScalarOpFunction(atan)
UnaryScalarOpFunction(sinh)
UnaryScalarOpFunction(cosh)
UnaryScalarOpFunction(tanh)
UnaryScalarOpFunction(asinh)
UnaryScalarOpFunction(acosh)
UnaryScalarOpFunction(atanh)
UnaryScalarOpFunction(erf)
UnaryScalarOpFunction(erfc)
UnaryScalarOpFunction(lgamma)
UnaryScalarOpFunction(j0)
UnaryScalarOpFunction(j1)
UnaryScalarOpFunction(y0)
UnaryScalarOpFunction(y1)

#undef UnaryScalarOpFunction

BinaryFunction(stabilise);
BinaryFunction(pow);
BinaryFunction(integerPow);
BinaryFunction(integerRoot);

/*---------------------------------------------------------------------------*\
                            Class atan2Op Declaration
\*---------------------------------------------------------------------------*/

struct atan2Op
{
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return "atan2(" + a + ',' + b + ')';
    }

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
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& n, const word& s)
    {
        return "jn(" + n + ',' + s + ')';
    }

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
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& n, const word& s)
    {
        return "yn(" + n + ',' + s + ')';
    }

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
