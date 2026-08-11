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
    Unary and binary tensor functions and operators used by the expression
    templating system

\*---------------------------------------------------------------------------*/

#ifndef tensorExpressions_H
#define tensorExpressions_H

#include "expressions.H"
#include "dimensionSet.H"

#include "diagTensor.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Unary tensor functions

#define UnaryAllTensorOpFunction(RetType, Func)                                \
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
    static RetType value(const diagTensor& x)                                  \
    {                                                                          \
        return Func(x);                                                        \
    }                                                                          \
                                                                               \
    static RetType value(const sphericalTensor& x)                             \
    {                                                                          \
        return Func(x);                                                        \
    }                                                                          \
                                                                               \
    static RetType value(const symmTensor& x)                                  \
    {                                                                          \
        return Func(x);                                                        \
    }                                                                          \
                                                                               \
    static RetType value(const tensor& x)                                      \
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

#define UnaryTensorInvOpFunction(Func)                                         \
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
    static diagTensor value(const diagTensor& x)                               \
    {                                                                          \
        return Func(x);                                                        \
    }                                                                          \
                                                                               \
    static sphericalTensor value(const sphericalTensor& x)                     \
    {                                                                          \
        return Func(x);                                                        \
    }                                                                          \
                                                                               \
    static symmTensor value(const symmTensor& x)                               \
    {                                                                          \
        return Func(x);                                                        \
    }                                                                          \
                                                                               \
    static tensor value(const tensor& x)                                       \
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

#define UnarySymmTensorOpFunction(RetType, Func)                               \
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
    static RetType value(const symmTensor& x)                                  \
    {                                                                          \
        return Func(x);                                                        \
    }                                                                          \
                                                                               \
    static RetType value(const tensor& x)                                      \
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

#define UnaryTensorDevOpFunction(Func)                                         \
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
    static symmTensor value(const symmTensor& x)                               \
    {                                                                          \
        return Func(x);                                                        \
    }                                                                          \
                                                                               \
    static tensor value(const tensor& x)                                       \
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

UnaryAllTensorOpFunction(scalar, tr)
UnaryAllTensorOpFunction(sphericalTensor, sph)
UnaryAllTensorOpFunction(scalar, det)

UnaryOpFunction(diagTensor, tensor, diag)

UnaryTensorInvOpFunction(inv)

UnarySymmTensorOpFunction(symmTensor, symm)
UnarySymmTensorOpFunction(symmTensor, twoSymm)
UnarySymmTensorOpFunction(vector, eigenValues)
UnarySymmTensorOpFunction(tensor, eigenVectors)
UnaryOpFunction(symmTensor, symmTensor, innerSqr)

UnaryTensorDevOpFunction(dev)
UnaryTensorDevOpFunction(dev2)
UnaryTensorDevOpFunction(cof)

UnaryOpFunction(tensor, tensor, skew)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Unary Hodge dual function

struct hdualOp
{
    static const signed char precedence = 2;

    static word name(const word& x)
    {
        return '*' + x;
    }

    static vector value(const symmTensor& x)
    {
        return *x;
    }

    static vector value(const tensor& x)
    {
        return *x;
    }

    static tensor value(const vector& x)
    {
        return *x;
    }

    static dimensionSet value(const dimensionSet& x)
    {
        return hdual(x);
    }
};

UnaryOperator(*, hdual)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef UnaryAllTensorOpFunction
#undef UnaryTensorInvOpFunction
#undef UnarySymmTensorOpFunction
#undef UnaryTensorDevOpFunction

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
