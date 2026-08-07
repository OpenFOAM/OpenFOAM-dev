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

\*---------------------------------------------------------------------------*/

#include "UPstream.H"

#define TEMPLATE template<class Type>
#include "FieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void component
(
    Field<typename Field<Type>::cmptType>& res,
    const UList<Type>& f,
    const direction d
);


template<class Type>
void T(Field<Type>& res, const UList<Type>& f);


template<class Type, direction r>
void pow
(
    Field<typename powProduct<Type, r>::type>& res,
    const Field<Type>& vf
);

template<class Type, direction r>
void pow
(
    Field<typename powProduct<Type, r>::type>& res,
    const SubField<Type>& vf
);


template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const Field<Type>& f,
    typename powProduct<Type, r>::type
      = pTraits<typename powProduct<Type, r>::type>::zero
);

template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const SubField<Type>& f,
    typename powProduct<Type, r>::type
      = pTraits<typename powProduct<Type, r>::type>::zero
);

template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const tmp<Field<Type>>& tf,
    typename powProduct<Type, r>::type
      = pTraits<typename powProduct<Type, r>::type>::zero
);


template<class Type>
void sqr
(
    Field<typename outerProduct<Type, Type>::type>& res,
    const UList<Type>& vf
);

template<class Type>
void magSqr(Field<scalar>& res, const UList<Type>& f);


template<class Type>
void mag(Field<scalar>& res, const UList<Type>& f);


template<class Type>
void cmptMax(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f);


template<class Type>
void cmptMin(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f);


template<class Type>
void cmptAv(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f);


template<class Type>
void cmptMag(Field<Type>& res, const UList<Type>& f);


BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(Type, Type, -, negate)

BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
BINARY_OPERATOR(Type, Type, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2>                                             \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Type1, Type2>::type>& res,                          \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
);                                                                             \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Type, Form>::type>& res,                            \
    const UList<Type>& f1,                                                     \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
);                                                                             \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Form, Type>::type>& res,                            \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const UList<Type>& f1                                                      \
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
