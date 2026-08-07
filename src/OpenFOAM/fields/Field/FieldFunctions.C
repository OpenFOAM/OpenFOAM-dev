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

#include "SubField.H"
#include "PstreamReduceOps.H"
#include "FieldReuseFunctions.H"

#define TEMPLATE template<class Type>
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * */

template<class Type>
void component
(
    Field<typename Field<Type>::cmptType>& res,
    const UList<Type>& f,
    const direction d
)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_F_FUNC_S
    (
        cmptType, res, =, Type, f, .component, const direction, d
    )
}


template<class Type>
void T(Field<Type>& res, const UList<Type>& f)
{
    TFOR_ALL_F_OP_F_FUNC(Type, res, =, Type, f, T)
}


template<class Type, direction r>
void pow
(
    Field<typename powProduct<Type, r>::type>& res,
    const Field<Type>& vf
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    TFOR_ALL_F_OP_FUNC_F_S
    (
        powProductType, res, =, pow, Type, vf, powProductType,
        pTraits<powProductType>::zero
    )
}

template<class Type, direction r>
void pow
(
    Field<typename powProduct<Type, r>::type>& res,
    const SubField<Type>& vf
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    TFOR_ALL_F_OP_FUNC_F_S
    (
        powProductType, res, =, pow, Type, vf, powProductType,
        pTraits<powProductType>::zero
    )
}


template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const Field<Type>& f,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<Field<powProductType>> tRes
    (
        new Field<powProductType>(f.size())
    );
    pow<Type, r>(tRes.ref(), f);
    return tRes;
}

template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const SubField<Type>& f,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<Field<powProductType>> tRes
    (
        new Field<powProductType>(f.size())
    );
    pow<Type, r>(tRes.ref(), f);
    return tRes;
}

template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const tmp<Field<Type>>& tf,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<Field<powProductType>> tRes = reuseTmp<powProductType, Type>::New(tf);
    pow<Type, r>(tRes.ref(), tf());
    tf.clear();
    return tRes;
}


template<class Type>
void sqr
(
    Field<typename outerProduct<Type, Type>::type>& res,
    const UList<Type>& vf
)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    TFOR_ALL_F_OP_FUNC_F(outerProductType, res, =, sqr, Type, vf)
}


template<class Type>
void magSqr(Field<scalar>& res, const UList<Type>& f)
{
    TFOR_ALL_F_OP_FUNC_F(scalar, res, =, magSqr, Type, f)
}


template<class Type>
void mag(Field<scalar>& res, const UList<Type>& f)
{
    TFOR_ALL_F_OP_FUNC_F(scalar, res, =, mag, Type, f)
}


template<class Type>
void cmptMax(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptMax, Type, f)
}


template<class Type>
void cmptMin(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptMin, Type, f)
}


template<class Type>
void cmptAv(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptAv, Type, f)
}


template<class Type>
void cmptMag(Field<Type>& res, const UList<Type>& f)
{
    TFOR_ALL_F_OP_FUNC_F(Type, res, =, cmptMag, Type, f)
}


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
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    TFOR_ALL_F_OP_F_OP_F(productType, res, =, Type1, f1, Op, Type2, f2)        \
}                                                                              \
                                                                               \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Type, Form>::type>& res,                            \
    const UList<Type>& f1,                                                     \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    TFOR_ALL_F_OP_F_OP_S                                                       \
        (productType, res, =,Type, f1, Op, Form, static_cast<const Form&>(vs)) \
}                                                                              \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Form, Type>::type>& res,                            \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const UList<Type>& f1                                                      \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    TFOR_ALL_F_OP_S_OP_F                                                       \
        (productType, res, =,Form,static_cast<const Form&>(vs), Op, Type, f1)  \
}

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
