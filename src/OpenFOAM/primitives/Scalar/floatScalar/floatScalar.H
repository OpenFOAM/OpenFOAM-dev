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

Typedef
    Foam::floatScalar

Description
    Float precision floating point scalar type.

SourceFiles
    floatScalar.C

\*---------------------------------------------------------------------------*/

#ifndef floatScalar_H
#define floatScalar_H

#include "doubleFloat.H"
#include "direction.H"
#include "word.H"

#include <limits>
using std::numeric_limits;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef float floatScalar;

// Largest and smallest scalar values allowed in certain parts of the code
static const floatScalar floatScalarVGreat
    = numeric_limits<floatScalar>::max()/10;
static const floatScalar floatScalarVSmall
    = numeric_limits<floatScalar>::min();

static const floatScalar floatScalarSmall
    = numeric_limits<floatScalar>::epsilon();
static const floatScalar floatScalarGreat
    = 1.0/floatScalarSmall;

static const floatScalar floatScalarRootVGreat = std::sqrt(floatScalarVGreat);
static const floatScalar floatScalarRootVSmall = std::sqrt(floatScalarVSmall);

static const floatScalar floatScalarRootGreat = std::sqrt(floatScalarGreat);
static const floatScalar floatScalarRootSmall = std::sqrt(floatScalarSmall);

static const floatScalar floatScalarNaN
    = numeric_limits<floatScalar>::signaling_NaN();

//- Read whole of buf as a scalar. Return true if successful.
inline bool readScalar(const char* buf, floatScalar& s)
{
    char* endPtr;
    s = strtof(buf, &endPtr);

    return (*endPtr == '\0');
}

#define Scalar floatScalar
#define ScalarVGreat floatScalarVGreat
#define ScalarVSmall floatScalarVSmall
#define ScalarRootVGreat floatScalarRootVGreat
#define ScalarRootVSmall floatScalarRootVSmall
#define ScalarNaN floatScalarNaN
#define readScalar readFloatScalar

inline Scalar mag(const Scalar s)
{
    return ::fabsf(s);
}


MAXMIN(Scalar, Scalar, Scalar)
MAXMIN(Scalar, Scalar, int32_t)
MAXMIN(Scalar, int32_t, Scalar)
MAXMIN(Scalar, Scalar, int64_t)
MAXMIN(Scalar, int64_t, Scalar)


#define POW(type1, type2)                         \
                                                  \
inline Scalar pow(const type1 s, const type2 e)   \
{                                                 \
    return ::powf(Scalar(s), Scalar(e));          \
}

POW(Scalar, Scalar)
POW(Scalar, double)
POW(int32_t, Scalar)
POW(int64_t, Scalar)

#undef POW


#define transFunc(func)            \
inline Scalar func(const Scalar s) \
{                                  \
    return ::func##f(s);           \
}

#define cFunc(func)                \
inline Scalar func(const Scalar s) \
{                                  \
    return ::func##f(s);           \
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
    return ::hypotf(x, y);
}

inline Scalar atan2(const Scalar y, const Scalar x)
{
    return ::atan2f(y, x);
}

inline Scalar jn(const int n, const Scalar s)
{
    return ::jnf(n, s);
}

inline Scalar yn(const int n, const Scalar s)
{
    return ::ynf(n, s);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef Scalar
#undef ScalarVGreat
#undef ScalarVSmall
#undef ScalarRootVGreat
#undef ScalarRootVSmall
#undef ScalarNaN
#undef readScalar
#undef transFunc
#undef cFunc

#endif

// ************************************************************************* //
