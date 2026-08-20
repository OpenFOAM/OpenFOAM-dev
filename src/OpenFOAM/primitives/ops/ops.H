/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

Global
    Foam::*Op

Description
    Operation classes containing functions and operators used by
    templated reduction and distribution functions.

\*---------------------------------------------------------------------------*/

#ifndef ops_H
#define ops_H

#include "expressionName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class nopOp Declaration
\*---------------------------------------------------------------------------*/

struct nopOp
{
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static word name(const word& x)
    {
        return x;
    }

    template<class X>
    static X value(const X& x)
    {
        return x;
    }

    template<class X>
    auto operator()(const X& x) const
    {
        return x;
    }
};

/*---------------------------------------------------------------------------*\
                          Class promoteOp Declaration
\*---------------------------------------------------------------------------*/

struct promoteOp
{
    static const signed char precedence = 2;

    static word name(const word& x)
    {
        return x;
    }

    template<class X>
    static X value(const X& x)
    {
        return x;
    }

    template<class X>
    auto operator()(const X& x) const
    {
        return x;
    }
};

/*---------------------------------------------------------------------------*\
                          Class negateOp Declaration
\*---------------------------------------------------------------------------*/

struct negateOp
{
    static const signed char precedence = 2;

    static word name(const word& x)
    {
        return '-' + x;
    }

    template<class X>
    static auto value(const X& x)
    {
        return -x;
    }

    template<class X>
    auto operator()(const X& x) const
    {
        return -x;
    }
};

/*---------------------------------------------------------------------------*\
                          Class notOp Declaration
\*---------------------------------------------------------------------------*/

struct notOp
{
    static const signed char precedence = 2;

    static word name(const word& x)
    {
        return '!' + x;
    }

    template<class X>
    static auto value(const X& x)
    {
        return !x;
    }

    template<class X>
    auto operator()(const X& x) const
    {
        return !x;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Unary functions

#define UnaryFunctionOp(Func)                                                  \
    struct Func##Op                                                            \
    {                                                                          \
        static const signed char precedence =                                  \
            expression::Operation<void>::precedence;                           \
                                                                               \
        static word name(const word& x)                                        \
        {                                                                      \
            return #Func"(" + x + ')';                                         \
        }                                                                      \
                                                                               \
        template<class Type>                                                   \
        static auto value(const Type x)                                        \
        {                                                                      \
            return Func(x);                                                    \
        }                                                                      \
                                                                               \
        template<class Type>                                                   \
        auto operator()(const Type x) const                                    \
        {                                                                      \
            return Func(x);                                                    \
        }                                                                      \
    };

UnaryFunctionOp(sqr)
UnaryFunctionOp(mag)
UnaryFunctionOp(magSqr)
UnaryFunctionOp(cmptMax)
UnaryFunctionOp(cmptMin)
UnaryFunctionOp(cmptAv)
UnaryFunctionOp(cmptMag)

/*---------------------------------------------------------------------------*\
                            Class addOp Declaration
\*---------------------------------------------------------------------------*/

struct addOp
{
    static const signed char precedence = 0;

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

    template<class A, class B>
    auto operator()(const A& a, const B& b) const
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

// Alias for addOp for readability of summation reduction operations
using sumOp = addOp;

/*---------------------------------------------------------------------------*\
                          Class subtractOp Declaration
\*---------------------------------------------------------------------------*/

struct subtractOp
{
    static const signed char precedence = 0;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return a + '-' + b;
    }

    template<class A, class B>
    auto operator()(const A& a, const B& b) const
    {
        return a - b;
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return a - b;
    }
};

/*---------------------------------------------------------------------------*\
                          Class multiplyOp Declaration
\*---------------------------------------------------------------------------*/

struct multiplyOp
{
    static const signed char precedence = 1;

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

    template<class A, class B>
    auto operator()(const A& a, const B& b) const
    {
        return a*b;
    }
};

/*---------------------------------------------------------------------------*\
                            Class divideOp Declaration
\*---------------------------------------------------------------------------*/

struct divideOp
{
    static const signed char precedence = 1;

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

    template<class A, class B>
    auto operator()(const A& a, const B& b) const
    {
        return a/b;
    }
};

/*---------------------------------------------------------------------------*\
                            Class innerOp Declaration
\*---------------------------------------------------------------------------*/

struct innerOp
{
    static const signed char precedence = 1;

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

/*---------------------------------------------------------------------------*\
                            Class crossOp Declaration
\*---------------------------------------------------------------------------*/

struct crossOp
{
    static const signed char precedence = 1;

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

/*---------------------------------------------------------------------------*\
                            Class stabiliseOp Declaration
\*---------------------------------------------------------------------------*/

struct stabiliseOp
{
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return "stabilise(" + a + ',' + b + ')';
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return stabilise(a, b);
    }

    scalar operator()(const scalar a, const scalar b) const
    {
        return stabilise(a, b);
    }
};

/*---------------------------------------------------------------------------*\
                           Class maxOp Declaration
\*---------------------------------------------------------------------------*/

struct maxOp
{
    static const signed char precedence =
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

    template<class A, class B>
    auto operator()(const A& a, const B& b) const
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

/*---------------------------------------------------------------------------*\
                             Class minOp Declaration
\*---------------------------------------------------------------------------*/

struct minOp
{
    static const signed char precedence =
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

    template<class A, class B>
    auto operator()(const A& a, const B& b) const
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

/*---------------------------------------------------------------------------*\
                             Class minModOp Declaration
\*---------------------------------------------------------------------------*/

struct minModOp
{
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return "minMod(" + a + ',' + b + ')';
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return minMod(a, b);
    }

    template<class A, class B>
    auto operator()(const A& a, const B& b) const
    {
        return minMod(a, b);
    }

    static word name(const word& x)
    {
        return "minMod(" + x + ')';
    }

    template<class X>
    static const X& identity()
    {
        return pTraits<X>::max;
    }
};

/*---------------------------------------------------------------------------*\
                          Class cmptMultiplyOp Declaration
\*---------------------------------------------------------------------------*/

struct cmptMultiplyOp
{
    static const signed char precedence = 1;

    static const bool associative = true;

    static word name(const word& a, const word& b)
    {
        return "cmptMultiply(" + a + ',' + b + ')';
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return cmptMultiply(a, b);
    }

    template<class A, class B>
    auto operator()(const A& a, const B& b) const
    {
        return cmptMultiply(a, b);
    }
};

/*---------------------------------------------------------------------------*\
                         Class cmptDivideOp Declaration
\*---------------------------------------------------------------------------*/

struct cmptDivideOp
{
    static const signed char precedence = 1;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return "cmptDivide(" + a + ',' + b + ')';
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return cmptDivide(a, b);
    }

    template<class A, class B>
    auto operator()(const A& a, const B& b) const
    {
        return cmptDivide(a, b);
    }
};

/*---------------------------------------------------------------------------*\
                           Class cmptPowOp Declaration
\*---------------------------------------------------------------------------*/

struct cmptPowOp
{
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& a, const word& b)
    {
        return "cmptPow(" + a + ',' + b + ')';
    }

    template<class A, class B>
    static auto value(const A& a, const B& b)
    {
        return cmptPow(a, b);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define EqOp(opName, op)                                                       \
    struct opName##Op                                                          \
    {                                                                          \
        template<class T1, class T2>                                           \
        void operator()(T1& x, const T2& y) const                              \
        {                                                                      \
            op;                                                                \
        }                                                                      \
    };

EqOp(eq, x = y)
EqOp(eqNeg, x = -y)
EqOp(eqSqr, x = sqr(y))
EqOp(eqMag, x = mag(y))

EqOp(addEq, x += y)
EqOp(subtractEq, x -= y)
EqOp(multiplyEq, x *= y)
EqOp(divideEq, x /= y)
EqOp(addEqMagSqr, x += magSqr(y))
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

#define Op(opName, op)                                                         \
    struct opName##Op                                                          \
    {                                                                          \
        template<class T1, class T2>                                           \
        auto operator()(const T1& x, const T2& y) const                        \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };

Op(minMagSqr, (magSqr(x) <= magSqr(y) ? x : y))
Op(maxMagSqr, (magSqr(x) >= magSqr(y) ? x : y))
Op(and, x && y)
Op(or, x || y)
Op(equal, x == y)
Op(less, x < y)
Op(lessEq, x <= y)
Op(greater, x > y)
Op(greaterEq, x >= y)

#undef Op

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
