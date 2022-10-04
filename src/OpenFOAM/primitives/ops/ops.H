/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

InClass
    Foam::Pstream

Description
    Combination-Reduction operation for a parallel run.

    The information from all nodes is collected on the master node,
    combined using the given combination function and the result is
    broadcast to all nodes

\*---------------------------------------------------------------------------*/

#ifndef ops_H
#define ops_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define EqOp(opName, op)                                                       \
                                                                               \
    template<class T1, class T2>                                               \
    class opName##Op2                                                          \
    {                                                                          \
    public:                                                                    \
                                                                               \
        void operator()(T1& x, const T2& y) const                              \
        {                                                                      \
            op;                                                                \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T>                                                          \
    class opName##Op                                                           \
    {                                                                          \
    public:                                                                    \
                                                                               \
        void operator()(T& x, const T& y) const                                \
        {                                                                      \
            op;                                                                \
        }                                                                      \
    };

EqOp(eq, x = y)
EqOp(eqNeg, x = -y)
EqOp(eqSqr, x = sqr(y))
EqOp(eqMag, x = mag(y))

EqOp(plusEq, x += y)
EqOp(minusEq, x -= y)
EqOp(multiplyEq, x *= y)
EqOp(divideEq, x /= y)
EqOp(plusEqMagSqr, x += magSqr(y))
EqOp(maxEq, x = max(x, y))
EqOp(minEq, x = min(x, y))
EqOp(minMagSqrEq, x = (magSqr(x) <= magSqr(y) ? x : y))
EqOp(maxMagSqrEq, x = (magSqr(x) >= magSqr(y) ? x : y))
EqOp(andEq, x = (x && y))
EqOp(orEq, x = (x || y))
EqOp(notEq, x = (x != y))

EqOp(nopEq, (void)x)

#undef EqOp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if __GNUC__
#define WARNRETURN __attribute__((warn_unused_result))
#else
#define WARNRETURN
#endif

#define Op(opName, op)                                                         \
                                                                               \
    template<class T, class T1, class T2>                                      \
    class opName##Op3                                                          \
    {                                                                          \
    public:                                                                    \
                                                                               \
        T operator()(const T1& x, const T2& y) const WARNRETURN                \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T1, class T2>                                               \
    class opName##Op2                                                          \
    {                                                                          \
    public:                                                                    \
                                                                               \
        T1 operator()(const T1& x, const T2& y) const WARNRETURN               \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T>                                                          \
    class opName##Op                                                           \
    {                                                                          \
    public:                                                                    \
                                                                               \
        T operator()(const T& x, const T& y) const WARNRETURN                  \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };

Op(sum, x + y)
Op(plus, x + y)
Op(minus, x - y)
Op(multiply, x * y)
Op(divide, x / y)
Op(cmptMultiply, cmptMultiply(x, y))
Op(cmptPow, cmptPow(x, y))
Op(cmptDivide, cmptDivide(x, y))
Op(stabilise, stabilise(x, y))
Op(max, max(x, y))
Op(min, min(x, y))
Op(minMagSqr, (magSqr(x) <= magSqr(y) ? x : y))
Op(maxMagSqr, (magSqr(x) >= magSqr(y) ? x : y))
Op(minMod, minMod(x, y))
Op(and, x && y)
Op(or, x || y)
Op(not, x != y)
Op(equal, x == y)
Op(less, x < y)
Op(lessEq, x <= y)
Op(greater, x > y)
Op(greaterEq, x >= y)

#undef Op
#undef WARNRETURN

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
