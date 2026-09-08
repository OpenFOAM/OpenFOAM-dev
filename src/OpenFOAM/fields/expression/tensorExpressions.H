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

#include "diagTensor.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "tensor.H"
#include "expressions.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Unary tensor functions

UnaryOpFunction(tr)
UnaryOpFunction(sph)
UnaryOpFunction(det)
UnaryOpFunction(diag)
UnaryOpFunction(T)
UnaryOpFunction(inv)
UnaryOpFunction(symm)
UnaryOpFunction(twoSymm)
UnaryOpFunction(eigenValues)
UnaryOpFunction(eigenVectors)
UnaryOpFunction(innerSqr)
UnaryOpFunction(dev)
UnaryOpFunction(dev2)
UnaryOpFunction(cof)
UnaryOpFunction(skew)

/*---------------------------------------------------------------------------*\
                            Class hdualOp Declaration
\*---------------------------------------------------------------------------*/

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

    template<class Type>
    static auto value(const Type& x)
    {
        return hdual(x);
    }
};

// Unary Hodge dual function
UnaryOperator(*, hdual)

/*---------------------------------------------------------------------------*\
                            Class ddotOp Declaration
\*---------------------------------------------------------------------------*/

struct ddotOp
{
    static const signed char precedence = 1;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return a + "&&" + b;
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return a && b;
    }
};

// Inner product (ddot) operator
BinaryOperator(&&, ddot);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
