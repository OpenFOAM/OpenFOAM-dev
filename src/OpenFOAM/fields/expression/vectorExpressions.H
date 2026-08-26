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
    Unary and binary vector functions and operators used by the expression
    templating system

\*---------------------------------------------------------------------------*/

#ifndef vectorExpressions_H
#define vectorExpressions_H

#include "expressions.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Unary vector functions

UnaryFunction(cmptMax)
UnaryFunction(cmptMin)
UnaryFunction(cmptAv)
UnaryFunction(cmptMag)

#define UnaryVectorOpFunction(func)                                            \
    UnaryOpFunction(vector, vector, func)

UnaryVectorOpFunction(normalised)
UnaryVectorOpFunction(perpendicular)

#undef UnaryVectorOpFunction

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Binary vector functions

BinaryFunction(cmptMultiply);
BinaryFunction(cmptDivide);

/*---------------------------------------------------------------------------*\
                         Class componentOp Declaration
\*---------------------------------------------------------------------------*/

struct componentOp
{
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& vs, const word& d)
    {
        return "component(" + vs + ',' + d + ')';
    }

    template<class VSType, class Dtype>
    static auto value(const VSType& vs, const Dtype& d)
    {
        return component(vs, d);
    }
};

BinaryFunction(component);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
