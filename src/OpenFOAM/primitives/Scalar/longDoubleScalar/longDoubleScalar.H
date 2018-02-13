/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

Typedef
    Foam::longDoubleScalar

Description
    Lang double precision floating point scalar type.

SourceFiles
    longDoubleScalar.C

\*---------------------------------------------------------------------------*/

#ifndef longDoubleScalar_H
#define longDoubleScalar_H

#include "doubleFloat.H"
#include "direction.H"
#include "word.H"

#include <limits>
using std::numeric_limits;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef long double longDoubleScalar;

// Largest and smallest scalar values allowed in certain parts of the code.
static const longDoubleScalar longDoubleScalarVGreat
    = numeric_limits<long double>::max()/10;
static const longDoubleScalar longDoubleScalarVSmall
    = numeric_limits<double>::min();

static const longDoubleScalar longDoubleScalarSmall
    = 1e3*numeric_limits<long double>::epsilon();
static const longDoubleScalar longDoubleScalarGreat
    = 1.0/longDoubleScalarSmall;

static const longDoubleScalar longDoubleScalarRootVGreat
    = ::sqrtl(longDoubleScalarVGreat);
static const longDoubleScalar longDoubleScalarRootVSmall
    = ::sqrtl(longDoubleScalarVSmall);

static const longDoubleScalar longDoubleScalarRootGreat
    = ::sqrtl(longDoubleScalarGreat);
static const longDoubleScalar longDoubleScalarRootSmall
    = ::sqrtl(longDoubleScalarSmall);

//- Read whole of buf as a scalar. Return true if successful.
inline bool readScalar(const char* buf, longDoubleScalar& s)
{
    char* endPtr;
    s = strtold(buf, &endPtr);

    return (*endPtr == '\0');
}

#define Scalar longDoubleScalar
#define ScalarVGreat longDoubleScalarVGreat
#define ScalarVSmall longDoubleScalarVSmall
#define ScalarRootVGreat longDoubleScalarRootVGreat
#define ScalarRootVSmall longDoubleScalarRootVSmall
#define readScalar readLongDoubleScalar

inline Scalar mag(const Scalar s)
{
    return ::fabsl(s);
}

#define MAXMINPOW(retType, type1, type2)          \
                                                  \
MAXMIN(retType, type1, type2)                     \
                                                  \
inline double pow(const type1 s, const type2 e)   \
{                                                 \
    return ::powl(Scalar(s), Scalar(e));          \
}

MAXMINPOW(Scalar, Scalar, Scalar)
MAXMINPOW(Scalar, Scalar, int)
MAXMINPOW(Scalar, int, Scalar)
MAXMINPOW(Scalar, Scalar, long)
MAXMINPOW(Scalar, long, Scalar)
MAXMINPOW(Scalar, Scalar, float)
MAXMINPOW(Scalar, float, Scalar)
MAXMINPOW(Scalar, Scalar, double)
MAXMINPOW(Scalar, double, Scalar)

#undef MAXMINPOW

#define transFunc(func)            \
inline Scalar func(const Scalar s) \
{                                  \
    return ::func(s);              \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Scalar hypot(const Scalar x, const Scalar y)
{
    return ::hypotl(x, y);
}

inline Scalar atan2(const Scalar y, const Scalar x)
{
    return ::atan2l(y, x);
}

inline Scalar jn(const int n, const Scalar s)
{
    return ::jnl(n, s);
}

inline Scalar yn(const int n, const Scalar s)
{
    return ::ynl(n, s);
}

#undef Scalar
#undef ScalarVGreat
#undef ScalarVSmall
#undef ScalarRootVGreat
#undef ScalarRootVSmall
#undef readScalar
#undef transFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
