/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

namespace Foam
{

/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<class Type>
inline void component
(
    pointPatchField<typename Field<Type>::cmptType>& sf,
    const pointPatchField<Type>& f,
    const direction d
)
{}


template<class Type>
inline void T
(
    pointPatchField<Type>& f1,
    const pointPatchField<Type>& f2
)
{}


template<class Type, direction r>
inline void pow
(
    Field<typename powProduct<Type, r>::type>& f,
    const pointPatchField<Type>& vf
)
{}


template<class Type>
inline void sqr
(
    Field<typename outerProduct<Type, Type>::type>& f,
    const pointPatchField<Type>& vf
)
{}


template<class Type>
inline void magSqr
(
    pointPatchField<scalar>& sf,
    const pointPatchField<Type>& f
)
{}


template<class Type>
inline void mag
(
    pointPatchField<scalar>& sf,
    const pointPatchField<Type>& f
)
{}


template<class Type>
inline void cmptAv
(
    pointPatchField<typename Field<Type>::cmptType>& cf,
    const pointPatchField<Type>& f
)
{}


template<class Type>
inline void cmptMag
(
    pointPatchField<Type>& cf,
    const pointPatchField<Type>& f
)
{}


#define BINARY_FUNCTION(func)                                                  \
                                                                               \
template<class Type>                                                           \
inline void func                                                               \
(                                                                              \
    pointPatchField<Type>& f,                                                  \
    const pointPatchField<Type>& f1,                                           \
    const pointPatchField<Type>& f2                                            \
)                                                                              \
{}                                                                             \
                                                                               \
template<class Type>                                                           \
inline void func                                                               \
(                                                                              \
    pointPatchField<Type>& f,                                                  \
    const pointPatchField<Type>& f1,                                           \
    const Type& s                                                              \
)                                                                              \
{}

BINARY_FUNCTION(max)
BINARY_FUNCTION(min)
BINARY_FUNCTION(cmptMultiply)
BINARY_FUNCTION(cmptDivide)


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

#define UNARY_OPERATOR(op, opFunc)                                             \
                                                                               \
template<class Type>                                                           \
inline void opFunc                                                             \
(                                                                              \
    pointPatchField<Type>& f,                                                  \
    const pointPatchField<Type>& f1                                            \
)                                                                              \
{}

UNARY_OPERATOR(-, negate)

#define BINARY_OPERATOR(Type1, Type2, op, opFunc)                              \
                                                                               \
template<class Type>                                                           \
inline void opFunc                                                             \
(                                                                              \
    pointPatchField<Type>& f,                                                  \
    const pointPatchField<Type1>& f1,                                          \
    const pointPatchField<Type2>& f2                                           \
)                                                                              \
{}

BINARY_OPERATOR(scalar, Type, *, multiply)
BINARY_OPERATOR(Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, /, divide)

#define BINARY_TYPE_OPERATOR_SF(TYPE, op, opFunc)                              \
                                                                               \
template<class Type>                                                           \
inline void opFunc                                                             \
(                                                                              \
    pointPatchField<Type>& f,                                                  \
    const TYPE& s,                                                             \
    const pointPatchField<Type>& f1                                            \
)                                                                              \
{}


#define BINARY_TYPE_OPERATOR_FS(TYPE, op, opFunc)                              \
                                                                               \
template<class Type>                                                           \
inline void opFunc                                                             \
(                                                                              \
    pointPatchField<Type>& f,                                                  \
    const pointPatchField<Type>& f1,                                           \
    const TYPE& s                                                              \
)                                                                              \
{}


BINARY_TYPE_OPERATOR_SF(scalar, *, multiply)
BINARY_TYPE_OPERATOR_FS(scalar, *, multiply)
BINARY_TYPE_OPERATOR_FS(scalar, /, divide)


#define PRODUCT_OPERATOR(product, op, opFunc)                                  \
                                                                               \
template                                                                       \
<                                                                              \
    class Type1,                                                               \
    class Type2                                                                \
>                                                                              \
inline void opFunc                                                             \
(                                                                              \
    pointPatchField                                                            \
    <typename product<Type1, Type2>::type>& f,                                 \
    const pointPatchField<Type1>& f1,                                          \
    const pointPatchField<Type2>& f2                                           \
)                                                                              \
{}                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
inline void opFunc                                                             \
(                                                                              \
    pointPatchField                                                            \
    <typename product<Type, Form>::type>& f,                                   \
    const pointPatchField<Type>& f1,                                           \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{}                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type                                                                 \
>                                                                              \
inline void opFunc                                                             \
(                                                                              \
    pointPatchField                                                            \
    <typename product<Form, Type>::type>& f,                                   \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const pointPatchField<Type>& f1                                            \
)                                                                              \
{}

PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


inline void hdual
(
    pointPatchField<vector>&,
    const pointPatchField<tensor>&
)
{}

inline void hdual
(
    pointPatchField<tensor>&,
    const pointPatchField<vector>&
)
{}

inline void diag
(
    pointPatchField<vector>&,
    const pointPatchField<tensor>&
)
{}

inline void tr
(
    pointPatchField<scalar>&,
    const pointPatchField<tensor>&
)
{}

inline void dev
(
    pointPatchField<tensor>&,
    const pointPatchField<tensor>&
)
{}

inline void dev2
(
    pointPatchField<tensor>&,
    const pointPatchField<tensor>&
)
{}

inline void det
(
    pointPatchField<scalar>&,
    const pointPatchField<tensor>&
)
{}

inline void inv
(
    pointPatchField<tensor>&,
    const pointPatchField<tensor>&
)
{}

inline void symm
(
    pointPatchField<tensor>&,
    const pointPatchField<tensor>&
)
{}

inline void twoSymm
(
    pointPatchField<tensor>&,
    const pointPatchField<tensor>&
)
{}

inline void skew
(
    pointPatchField<tensor>&,
    const pointPatchField<tensor>&
)
{}

inline void eigenValues
(
    pointPatchField<vector>&,
    const pointPatchField<tensor>&
)
{}

inline void eigenVectors
(
    pointPatchField<tensor>&,
    const pointPatchField<tensor>&
)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
