/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "scalarFieldField.H"

#define TEMPLATE template<template<class> class Field, class Type>
#include "FieldFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<template<class> class Field, class Type>
void component
(
    FieldField<Field, typename FieldField<Field, Type>::cmptType>& sf,
    const FieldField<Field, Type>& f,
    const direction d
);

template<template<class> class Field, class Type>
void T(FieldField<Field, Type>& f1, const FieldField<Field, Type>& f2);


template<template<class> class Field, class Type, int r>
void pow
(
    FieldField<Field, typename powProduct<Type, r>::type>& f,
    const FieldField<Field, Type>& vf
);

template<template<class> class Field, class Type, int r>
tmp<FieldField<Field, typename powProduct<Type, r>::type> >
pow
(
    const FieldField<Field, Type>& f, typename powProduct<Type, r>::type
        = pTraits<typename powProduct<Type, r>::type>::zero
);

template<template<class> class Field, class Type, int r>
tmp<FieldField<Field, typename powProduct<Type, r>::type> >
pow
(
    const tmp<FieldField<Field, Type> >& tf, typename powProduct<Type, r>::type
        = pTraits<typename powProduct<Type, r>::type>::zero
);


template<template<class> class Field, class Type>
void sqr
(
    FieldField<Field, typename outerProduct<Type, Type>::type>& f,
    const FieldField<Field, Type>& vf
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename outerProduct<Type, Type>::type> >
sqr(const FieldField<Field, Type>& f);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename outerProduct<Type, Type>::type> >
sqr(const tmp<FieldField<Field, Type> >& tf);

template<template<class> class Field, class Type>
void magSqr(FieldField<Field, scalar>& sf, const FieldField<Field, Type>& f);

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > magSqr(const FieldField<Field, Type>& f);

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > magSqr(const tmp<FieldField<Field, Type> >& tf);

template<template<class> class Field, class Type>
void mag(FieldField<Field, scalar>& sf, const FieldField<Field, Type>& f);

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > mag(const FieldField<Field, Type>& f);

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > mag(const tmp<FieldField<Field, Type> >& tf);


template<template<class> class Field, class Type>
void cmptMax
(
    FieldField<Field, typename FieldField<Field, Type>::cmptType>& cf,
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename FieldField<Field, Type>::cmptType> > cmptMax
(
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename FieldField<Field, Type>::cmptType> > cmptMax
(
    const tmp<FieldField<Field, Type> >& tf
);


template<template<class> class Field, class Type>
void cmptMin
(
    FieldField<Field, typename FieldField<Field, Type>::cmptType>& cf,
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename FieldField<Field, Type>::cmptType> > cmptMin
(
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename FieldField<Field, Type>::cmptType> > cmptMin
(
    const tmp<FieldField<Field, Type> >& tf
);


template<template<class> class Field, class Type>
void cmptAv
(
    FieldField<Field, typename FieldField<Field, Type>::cmptType>& cf,
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename FieldField<Field, Type>::cmptType> > cmptAv
(
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename FieldField<Field, Type>::cmptType> > cmptAv
(
    const tmp<FieldField<Field, Type> >& tf
);


template<template<class> class Field, class Type>
void cmptMag
(
    FieldField<Field, Type>& cf,
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type> > cmptMag
(
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type> > cmptMag
(
    const tmp<FieldField<Field, Type> >& tf
);


#define TMP_UNARY_FUNCTION(returnType, func)                                  \
                                                                              \
template<template<class> class Field, class Type>                             \
returnType func(const tmp<FieldField<Field, Type> >& tf1);

template<template<class> class Field, class Type>
Type max(const FieldField<Field, Type>& f);

TMP_UNARY_FUNCTION(Type, max)

template<template<class> class Field, class Type>
Type min(const FieldField<Field, Type>& f);

TMP_UNARY_FUNCTION(Type, min)

template<template<class> class Field, class Type>
Type sum(const FieldField<Field, Type>& f);

TMP_UNARY_FUNCTION(Type, sum)

template<template<class> class Field, class Type>
scalar sumMag(const FieldField<Field, Type>& f);

TMP_UNARY_FUNCTION(scalar, sumMag)

template<template<class> class Field, class Type>
Type average(const FieldField<Field, Type>& f);

TMP_UNARY_FUNCTION(Type, average)


#define G_UNARY_FUNCTION(returnType, gFunc, func, rFunc)                      \
                                                                              \
template<template<class> class Field, class Type>                             \
returnType gFunc(const FieldField<Field, Type>& f);                           \
TMP_UNARY_FUNCTION(returnType, gFunc)

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(scalar, gSumMag, sumMag, sum)

#undef G_UNARY_FUNCTION


template<template<class> class Field, class Type>
Type gAverage(const FieldField<Field, Type>& f);

TMP_UNARY_FUNCTION(Type, gAverage)

#undef TMP_UNARY_FUNCTION


BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

UNARY_OPERATOR(Type, Type, -, negate)

#ifndef __INTEL_COMPILER
BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
#endif
BINARY_OPERATOR(Type, Type, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, op, opFunc)                                  \
                                                                               \
template                                                                       \
<                                                                              \
    template<class> class Field1,                                              \
    template<class> class Field2,                                              \
    class Type1,                                                               \
    class Type2                                                                \
>                                                                              \
void opFunc                                                                    \
(                                                                              \
    FieldField<Field1, typename product<Type1, Type2>::type>& f,               \
    const FieldField<Field1, Type1>& f1,                                       \
    const FieldField<Field2, Type2>& f2                                        \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    template<class> class Field1,                                              \
    template<class> class Field2,                                              \
    class Type1,                                                               \
    class Type2                                                                \
>                                                                              \
tmp<FieldField<Field1, typename product<Type1, Type2>::type> >                 \
operator op                                                                    \
(                                                                              \
    const FieldField<Field1, Type1>& f1,                                       \
    const FieldField<Field2, Type2>& f2                                        \
);                                                                             \
                                                                               \
template<template<class> class Field, class Type1, class Type2>                \
tmp<FieldField<Field, typename product<Type1, Type2>::type> >                  \
operator op                                                                    \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2> >& tf2                                  \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    template<class> class Field1,                                              \
    template<class> class Field2,                                              \
    class Type1,                                                               \
    class Type2                                                                \
>                                                                              \
tmp<FieldField<Field1, typename product<Type1, Type2>::type> >                 \
operator op                                                                    \
(                                                                              \
    const FieldField<Field1, Type1>& f1,                                       \
    const tmp<FieldField<Field2, Type2> >& tf2                                 \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    template<class> class Field1,                                              \
    template<class> class Field2,                                              \
    class Type1,                                                               \
    class Type2                                                                \
>                                                                              \
tmp<FieldField<Field1, typename product<Type1, Type2>::type> >                 \
operator op                                                                    \
(                                                                              \
    const tmp<FieldField<Field1, Type1> >& tf1,                                \
    const FieldField<Field2, Type2>& f2                                        \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    template<class> class Field1,                                              \
    template<class> class Field2,                                              \
    class Type1,                                                               \
    class Type2                                                                \
>                                                                              \
tmp<FieldField<Field1, typename product<Type1, Type2>::type> >                 \
operator op                                                                    \
(                                                                              \
    const tmp<FieldField<Field1, Type1> >& tf1,                                \
    const tmp<FieldField<Field2, Type2> >& tf2                                 \
);                                                                             \
                                                                               \
template                                                                       \
<template<class> class Field, class Type, class Form, class Cmpt, int nCmpt>   \
void opFunc                                                                    \
(                                                                              \
    FieldField<Field, typename product<Type, Form>::type>& f,                  \
    const FieldField<Field, Type>& f1,                                         \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
);                                                                             \
                                                                               \
template                                                                       \
<template<class> class Field, class Type, class Form, class Cmpt, int nCmpt>   \
tmp<FieldField<Field, typename product<Type, Form>::type> >                    \
operator op                                                                    \
(                                                                              \
    const FieldField<Field, Type>& f1,                                         \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
);                                                                             \
                                                                               \
template                                                                       \
<template<class> class Field, class Type, class Form, class Cmpt, int nCmpt>   \
tmp<FieldField<Field, typename product<Type, Form>::type> >                    \
operator op                                                                    \
(                                                                              \
    const tmp<FieldField<Field, Type> >& tf1,                                  \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
);                                                                             \
                                                                               \
template                                                                       \
<template<class> class Field, class Form, class Cmpt, int nCmpt, class Type>   \
void opFunc                                                                    \
(                                                                              \
    FieldField<Field, typename product<Form, Type>::type>& f,                  \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const FieldField<Field, Type>& f1                                          \
);                                                                             \
                                                                               \
template                                                                       \
<template<class> class Field, class Form, class Cmpt, int nCmpt, class Type>   \
tmp<FieldField<Field, typename product<Form, Type>::type> >                    \
operator op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const FieldField<Field, Type>& f1                                          \
);                                                                             \
                                                                               \
template                                                                       \
<template<class> class Field, class Form, class Cmpt, int nCmpt, class Type>   \
tmp<FieldField<Field, typename product<Form, Type>::type> >                    \
operator op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<FieldField<Field, Type> >& tf1                                   \
);

PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
