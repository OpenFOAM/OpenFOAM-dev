/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Specialisation of FieldField\<T\> for scalar.

\*---------------------------------------------------------------------------*/

#include "scalarFieldField.H"

#define TEMPLATE template<template<class> class Field>
#include "FieldFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class Field>
void stabilise
(
    FieldField<Field, scalar>& f,
    const FieldField<Field, scalar>& f1,
    const scalar s
)
{
    forAll(f, i)
    {
        stabilise(f[i], f1[i], s);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > stabilise
(
    const FieldField<Field, scalar>& f1,
    const scalar s
)
{
    tmp<FieldField<Field, scalar> > tf
    (
        FieldField<Field, scalar>::NewCalculatedType(f1)
    );
    stabilise(tf(), f1, s);
    return tf;
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > stabilise
(
    const tmp<FieldField<Field, scalar> >& tf1,
    const scalar s
)
{
    tmp<FieldField<Field, scalar> > tf(tf1.ptr());
    stabilise(tf(), tf(), s);
    return tf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(scalar, scalar, scalar, +, add)
BINARY_TYPE_OPERATOR(scalar, scalar, scalar, -, subtract)

BINARY_OPERATOR(scalar, scalar, scalar, *, multiply)
BINARY_OPERATOR(scalar, scalar, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(scalar, scalar, scalar, /, divide)

BINARY_FUNCTION(scalar, scalar, scalar, pow)
BINARY_TYPE_FUNCTION(scalar, scalar, scalar, pow)

BINARY_FUNCTION(scalar, scalar, scalar, atan2)
BINARY_TYPE_FUNCTION(scalar, scalar, scalar, atan2)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, scalar, pow3)
UNARY_FUNCTION(scalar, scalar, pow4)
UNARY_FUNCTION(scalar, scalar, pow5)
UNARY_FUNCTION(scalar, scalar, pow6)
UNARY_FUNCTION(scalar, scalar, pow025)
UNARY_FUNCTION(scalar, scalar, sqrt)
UNARY_FUNCTION(scalar, scalar, cbrt)
UNARY_FUNCTION(scalar, scalar, sign)
UNARY_FUNCTION(scalar, scalar, pos)
UNARY_FUNCTION(scalar, scalar, neg)
UNARY_FUNCTION(scalar, scalar, exp)
UNARY_FUNCTION(scalar, scalar, log)
UNARY_FUNCTION(scalar, scalar, log10)
UNARY_FUNCTION(scalar, scalar, sin)
UNARY_FUNCTION(scalar, scalar, cos)
UNARY_FUNCTION(scalar, scalar, tan)
UNARY_FUNCTION(scalar, scalar, asin)
UNARY_FUNCTION(scalar, scalar, acos)
UNARY_FUNCTION(scalar, scalar, atan)
UNARY_FUNCTION(scalar, scalar, sinh)
UNARY_FUNCTION(scalar, scalar, cosh)
UNARY_FUNCTION(scalar, scalar, tanh)
UNARY_FUNCTION(scalar, scalar, asinh)
UNARY_FUNCTION(scalar, scalar, acosh)
UNARY_FUNCTION(scalar, scalar, atanh)
UNARY_FUNCTION(scalar, scalar, erf)
UNARY_FUNCTION(scalar, scalar, erfc)
UNARY_FUNCTION(scalar, scalar, lgamma)
UNARY_FUNCTION(scalar, scalar, j0)
UNARY_FUNCTION(scalar, scalar, j1)
UNARY_FUNCTION(scalar, scalar, y0)
UNARY_FUNCTION(scalar, scalar, y1)


#define BesselFunc(func)                                                      \
                                                                              \
template<template<class> class Field>                                         \
void func                                                                     \
(                                                                             \
    FieldField<Field, scalar>& res,                                           \
    const int n,                                                              \
    const FieldField<Field, scalar>& sf                                       \
)                                                                             \
{                                                                             \
    forAll(res, i)                                                            \
    {                                                                         \
        func(res[i], n, sf[i]);                                               \
    }                                                                         \
}                                                                             \
                                                                              \
template<template<class> class Field>                                         \
tmp<FieldField<Field, scalar> > func                                          \
(                                                                             \
    const int n,                                                              \
    const FieldField<Field, scalar>& sf                                       \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, scalar> > tRes                                      \
    (                                                                         \
        FieldField<Field, scalar>::NewCalculatedType(sf)                      \
    );                                                                        \
    func(tRes(), n, sf);                                                      \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<template<class> class Field>                                         \
tmp<FieldField<Field, scalar> > func                                          \
(                                                                             \
    const int n,                                                              \
    const tmp<FieldField<Field, scalar> >& tsf                                \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, scalar> > tRes                                      \
    (                                                                         \
        reuseTmpFieldField<Field, scalar, scalar>::New(tsf)                   \
    );                                                                        \
    func(tRes(), n, tsf());                                                   \
    reuseTmpFieldField<Field, scalar, scalar>::clear(tsf);                    \
    return tRes;                                                              \
}

BesselFunc(jn)
BesselFunc(yn)

#undef BesselFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
