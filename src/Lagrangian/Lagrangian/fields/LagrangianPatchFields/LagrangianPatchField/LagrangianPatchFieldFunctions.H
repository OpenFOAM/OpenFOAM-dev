/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianPatchField.H"

#define FIELD_TYPE LagrangianPatchField

#define COMMA ,

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
inline void component
(
    FIELD_TYPE<typename FIELD_TYPE<Type>::cmptType>& sf,
    const FIELD_TYPE<Type>& f,
    const direction d
)
{}


#define UNARY_FUNCTION(Template, Type, Type1, func)                            \
                                                                               \
Template                                                                       \
inline void func                                                               \
(                                                                              \
    FIELD_TYPE<Type>& f,                                                       \
    const FIELD_TYPE<Type1>& f1                                                \
)                                                                              \
{}

UNARY_FUNCTION(template<class Type>, scalar, Type, magSqr)
UNARY_FUNCTION(template<class Type>, scalar, Type, mag)
UNARY_FUNCTION
(
    template<class Type>,
    typename FIELD_TYPE<Type>::cmptType,
    Type,
    cmptAv
)
UNARY_FUNCTION(template<class Type>, Type, Type, cmptMag)
UNARY_FUNCTION
(
    template<class Type COMMA direction r>,
    typename powProduct<Type COMMA r>::type,
    Type,
    pow
)
UNARY_FUNCTION
(
    template<class Type>,
    typename outerProduct<Type COMMA Type>::type,
    Type,
    sqr
)
UNARY_FUNCTION(, scalar, scalar, pow3)
UNARY_FUNCTION(, scalar, scalar, pow4)
UNARY_FUNCTION(, scalar, scalar, pow5)
UNARY_FUNCTION(, scalar, scalar, pow6)
UNARY_FUNCTION(, scalar, scalar, pow025)
UNARY_FUNCTION(, scalar, scalar, sqrt)
UNARY_FUNCTION(, scalar, scalar, cbrt)
UNARY_FUNCTION(template<class Type>, Type, Type, T)
UNARY_FUNCTION(, vector, tensor, hdual)
UNARY_FUNCTION(, tensor, vector, hdual)
UNARY_FUNCTION(, vector, tensor, diag)
UNARY_FUNCTION(, scalar, tensor, tr)
UNARY_FUNCTION(, tensor, tensor, dev)
UNARY_FUNCTION(, tensor, tensor, dev2)
UNARY_FUNCTION(, scalar, tensor, det)
UNARY_FUNCTION(, tensor, tensor, inv)
UNARY_FUNCTION(, tensor, tensor, symm)
UNARY_FUNCTION(, tensor, tensor, twoSymm)
UNARY_FUNCTION(, tensor, tensor, skew)
UNARY_FUNCTION(, vector, tensor, eigenValues)
UNARY_FUNCTION(, tensor, tensor, eigenVectors)

#undef UNARY_FUNCTION


#define BINARY_FUNCTION(func)                                                  \
                                                                               \
template<class Type>                                                           \
inline void func                                                               \
(                                                                              \
    FIELD_TYPE<Type>& f,                                                       \
    const FIELD_TYPE<Type>& f1,                                                \
    const FIELD_TYPE<Type>& f2                                                 \
)                                                                              \
{}                                                                             \
                                                                               \
template<class Type>                                                           \
inline void func                                                               \
(                                                                              \
    FIELD_TYPE<Type>& f,                                                       \
    const FIELD_TYPE<Type>& f1,                                                \
    const Type& s                                                              \
)                                                                              \
{}

BINARY_FUNCTION(max)
BINARY_FUNCTION(min)
BINARY_FUNCTION(cmptMultiply)
BINARY_FUNCTION(cmptDivide)

#undef BINARY_FUNCTION


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

#define UNARY_OPERATOR(op, opFunc)                                             \
                                                                               \
template<class Type>                                                           \
inline void opFunc(FIELD_TYPE<Type>& f, const FIELD_TYPE<Type>& f1)            \
{}

UNARY_OPERATOR(-, negate)

#undef UNARY_OPERATOR


#define BINARY_OPERATOR(Template, Type, Type1, Type2, op, opFunc)              \
                                                                               \
Template                                                                       \
inline void opFunc                                                             \
(                                                                              \
    FIELD_TYPE<Type>& f,                                                       \
    const FIELD_TYPE<Type1>& f1,                                               \
    const FIELD_TYPE<Type2>& f2                                                \
)                                                                              \
{}

BINARY_OPERATOR(template<class Type>, Type, scalar, Type, *, multiply)
BINARY_OPERATOR(template<class Type>, Type, Type, scalar, *, multiply)
BINARY_OPERATOR(, scalar, scalar, scalar, *, multiply)
BINARY_OPERATOR(template<class Type>, Type, Type, scalar, /, divide)

#undef BINARY_OPERATOR


#define BINARY_TYPE_OPERATOR_SF(TYPE, op, opFunc)                              \
                                                                               \
template<class Type>                                                           \
inline void opFunc                                                             \
(                                                                              \
    FIELD_TYPE<Type>& f,                                                       \
    const TYPE& s,                                                             \
    const FIELD_TYPE<Type>& f1                                                 \
)                                                                              \
{}

#define BINARY_TYPE_OPERATOR_FS(TYPE, op, opFunc)                              \
                                                                               \
template<class Type>                                                           \
inline void opFunc                                                             \
(                                                                              \
    FIELD_TYPE<Type>& f,                                                       \
    const FIELD_TYPE<Type>& f1,                                                \
    const TYPE& s                                                              \
)                                                                              \
{}

BINARY_TYPE_OPERATOR_SF(scalar, *, multiply)
BINARY_TYPE_OPERATOR_FS(scalar, *, multiply)
BINARY_TYPE_OPERATOR_FS(scalar, /, divide)

#undef BINARY_TYPE_OPERATOR_SF
#undef BINARY_TYPE_OPERATOR_FS


#define PRODUCT_OPERATOR(product, op, opFunc)                                  \
                                                                               \
template<class Type1, class Type2>                                             \
inline void opFunc                                                             \
(                                                                              \
    FIELD_TYPE<typename product<Type1, Type2>::type>& f,                       \
    const FIELD_TYPE<Type1>& f1,                                               \
    const FIELD_TYPE<Type2>& f2                                                \
)                                                                              \
{}                                                                             \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
inline void opFunc                                                             \
(                                                                              \
    FIELD_TYPE<typename product<Type, Form>::type>& f,                         \
    const FIELD_TYPE<Type>& f1,                                                \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{}                                                                             \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
inline void opFunc                                                             \
(                                                                              \
    FIELD_TYPE<typename product<Form, Type>::type>& f,                         \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const FIELD_TYPE<Type>& f1                                                 \
)                                                                              \
{}

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

#undef COMMA

#undef FIELD_TYPE

// ************************************************************************* //
