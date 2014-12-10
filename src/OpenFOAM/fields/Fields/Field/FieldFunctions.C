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


template<class Type, int r>
void pow
(
    Field<typename powProduct<Type, r>::type>& res,
    const UList<Type>& vf
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    TFOR_ALL_F_OP_FUNC_F_S
    (
        powProductType, res, =, pow, Type, vf, powProductType,
        pTraits<powProductType>::zero
    )
}

template<class Type, int r>
tmp<Field<typename powProduct<Type, r>::type> >
pow
(
    const UList<Type>& f,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<Field<powProductType> > tRes
    (
        new Field<powProductType>(f.size())
    );
    pow<Type, r>(tRes(), f);
    return tRes;
}

template<class Type, int r>
tmp<Field<typename powProduct<Type, r>::type> >
pow
(
    const tmp<Field<Type> >& tf,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<Field<powProductType> > tRes = reuseTmp<powProductType, Type>::New(tf);
    pow<Type, r>(tRes(), tf());
    reuseTmp<powProductType, Type>::clear(tf);
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
tmp<Field<typename outerProduct<Type, Type>::type> >
sqr(const UList<Type>& f)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<Field<outerProductType> > tRes
    (
        new Field<outerProductType>(f.size())
    );
    sqr(tRes(), f);
    return tRes;
}

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type> >
sqr(const tmp<Field<Type> >& tf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<Field<outerProductType> > tRes =
        reuseTmp<outerProductType, Type>::New(tf);
    sqr(tRes(), tf());
    reuseTmp<outerProductType, Type>::clear(tf);
    return tRes;
}


template<class Type>
void magSqr(Field<scalar>& res, const UList<Type>& f)
{
    TFOR_ALL_F_OP_FUNC_F(scalar, res, =, magSqr, Type, f)
}

template<class Type>
tmp<Field<scalar> > magSqr(const UList<Type>& f)
{
    tmp<Field<scalar> > tRes(new Field<scalar>(f.size()));
    magSqr(tRes(), f);
    return tRes;
}

template<class Type>
tmp<Field<scalar> > magSqr(const tmp<Field<Type> >& tf)
{
    tmp<Field<scalar> > tRes = reuseTmp<scalar, Type>::New(tf);
    magSqr(tRes(), tf());
    reuseTmp<scalar, Type>::clear(tf);
    return tRes;
}


template<class Type>
void mag(Field<scalar>& res, const UList<Type>& f)
{
    TFOR_ALL_F_OP_FUNC_F(scalar, res, =, mag, Type, f)
}

template<class Type>
tmp<Field<scalar> > mag(const UList<Type>& f)
{
    tmp<Field<scalar> > tRes(new Field<scalar>(f.size()));
    mag(tRes(), f);
    return tRes;
}

template<class Type>
tmp<Field<scalar> > mag(const tmp<Field<Type> >& tf)
{
    tmp<Field<scalar> > tRes = reuseTmp<scalar, Type>::New(tf);
    mag(tRes(), tf());
    reuseTmp<scalar, Type>::clear(tf);
    return tRes;
}


template<class Type>
void cmptMax(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptMax, Type, f)
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptMax(const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    tmp<Field<cmptType> > tRes(new Field<cmptType>(f.size()));
    cmptMax(tRes(), f);
    return tRes;
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptMax(const tmp<Field<Type> >& tf)
{
    typedef typename Field<Type>::cmptType cmptType;
    tmp<Field<cmptType> > tRes = reuseTmp<cmptType, Type>::New(tf);
    cmptMax(tRes(), tf());
    reuseTmp<cmptType, Type>::clear(tf);
    return tRes;
}


template<class Type>
void cmptMin(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptMin, Type, f)
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptMin(const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    tmp<Field<cmptType> > tRes(new Field<cmptType>(f.size()));
    cmptMin(tRes(), f);
    return tRes;
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptMin(const tmp<Field<Type> >& tf)
{
    typedef typename Field<Type>::cmptType cmptType;
    tmp<Field<cmptType> > tRes = reuseTmp<cmptType, Type>::New(tf);
    cmptMin(tRes(), tf());
    reuseTmp<cmptType, Type>::clear(tf);
    return tRes;
}


template<class Type>
void cmptAv(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptAv, Type, f)
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptAv(const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    tmp<Field<cmptType> > tRes(new Field<cmptType>(f.size()));
    cmptAv(tRes(), f);
    return tRes;
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptAv(const tmp<Field<Type> >& tf)
{
    typedef typename Field<Type>::cmptType cmptType;
    tmp<Field<cmptType> > tRes = reuseTmp<cmptType, Type>::New(tf);
    cmptAv(tRes(), tf());
    reuseTmp<cmptType, Type>::clear(tf);
    return tRes;
}


template<class Type>
void cmptMag(Field<Type>& res, const UList<Type>& f)
{
    TFOR_ALL_F_OP_FUNC_F(Type, res, =, cmptMag, Type, f)
}

template<class Type>
tmp<Field<Type> > cmptMag(const UList<Type>& f)
{
    tmp<Field<Type> > tRes(new Field<Type>(f.size()));
    cmptMag(tRes(), f);
    return tRes;
}

template<class Type>
tmp<Field<Type> > cmptMag(const tmp<Field<Type> >& tf)
{
    tmp<Field<Type> > tRes = reuseTmp<Type, Type>::New(tf);
    cmptMag(tRes(), tf());
    reuseTmp<Type, Type>::clear(tf);
    return tRes;
}


#define TMP_UNARY_FUNCTION(ReturnType, Func)                                  \
                                                                              \
template<class Type>                                                          \
ReturnType Func(const tmp<Field<Type> >& tf1)                                 \
{                                                                             \
    ReturnType res = Func(tf1());                                             \
    tf1.clear();                                                              \
    return res;                                                               \
}

template<class Type>
Type max(const UList<Type>& f)
{
    if (f.size())
    {
        Type Max(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S(Type, Max, =, max, Type, f, Type, Max)
        return Max;
    }
    else
    {
        return pTraits<Type>::min;
    }
}

TMP_UNARY_FUNCTION(Type, max)

template<class Type>
Type min(const UList<Type>& f)
{
    if (f.size())
    {
        Type Min(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S(Type, Min, =, min, Type, f, Type, Min)
        return Min;
    }
    else
    {
        return pTraits<Type>::max;
    }
}

TMP_UNARY_FUNCTION(Type, min)

template<class Type>
Type sum(const UList<Type>& f)
{
    if (f.size())
    {
        Type Sum = pTraits<Type>::zero;
        TFOR_ALL_S_OP_F(Type, Sum, +=, Type, f)
        return Sum;
    }
    else
    {
        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, sum)

template<class Type>
Type maxMagSqr(const UList<Type>& f)
{
    if (f.size())
    {
        Type Max(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S
        (
            Type,
            Max,
            =,
            maxMagSqrOp<Type>(),
            Type,
            f,
            Type,
            Max
        )
        return Max;
    }
    else
    {
        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, maxMagSqr)

template<class Type>
Type minMagSqr(const UList<Type>& f)
{
    if (f.size())
    {
        Type Min(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S
        (
            Type,
            Min,
            =,
            minMagSqrOp<Type>(),
            Type,
            f,
            Type,
            Min
        )
        return Min;
    }
    else
    {
        return pTraits<Type>::rootMax;
    }
}

TMP_UNARY_FUNCTION(Type, minMagSqr)

template<class Type>
scalar sumProd(const UList<Type>& f1, const UList<Type>& f2)
{
    if (f1.size() && (f1.size() == f2.size()))
    {
        scalar SumProd = 0.0;
        TFOR_ALL_S_OP_F_OP_F(scalar, SumProd, +=, Type, f1, &&, Type, f2)
        return SumProd;
    }
    else
    {
        return 0.0;
    }
}


template<class Type>
Type sumCmptProd(const UList<Type>& f1, const UList<Type>& f2)
{
    if (f1.size() && (f1.size() == f2.size()))
    {
        Type SumProd = pTraits<Type>::zero;
        TFOR_ALL_S_OP_FUNC_F_F
        (
            Type,
            SumProd,
            +=,
            cmptMultiply,
            Type,
            f1,
            Type,
            f2
        )
        return SumProd;
    }
    else
    {
        return pTraits<Type>::zero;
    }
}


template<class Type>
scalar sumSqr(const UList<Type>& f)
{
    if (f.size())
    {
        scalar SumSqr = 0.0;
        TFOR_ALL_S_OP_FUNC_F(scalar, SumSqr, +=, sqr, Type, f)
        return SumSqr;
    }
    else
    {
        return 0.0;
    }
}

TMP_UNARY_FUNCTION(scalar, sumSqr)

template<class Type>
scalar sumMag(const UList<Type>& f)
{
    if (f.size())
    {
        scalar SumMag = 0.0;
        TFOR_ALL_S_OP_FUNC_F(scalar, SumMag, +=, mag, Type, f)
        return SumMag;
    }
    else
    {
        return 0.0;
    }
}

TMP_UNARY_FUNCTION(scalar, sumMag)


template<class Type>
Type sumCmptMag(const UList<Type>& f)
{
    if (f.size())
    {
        Type SumMag = pTraits<Type>::zero;
        TFOR_ALL_S_OP_FUNC_F(scalar, SumMag, +=, cmptMag, Type, f)
        return SumMag;
    }
    else
    {
        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, sumCmptMag)

template<class Type>
Type average(const UList<Type>& f)
{
    if (f.size())
    {
        Type avrg = sum(f)/f.size();

        return avrg;
    }
    else
    {
        WarningIn("average(const UList<Type>&)")
            << "empty field, returning zero" << endl;

        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, average)


#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                      \
                                                                              \
template<class Type>                                                          \
ReturnType gFunc(const UList<Type>& f, const int comm)                        \
{                                                                             \
    ReturnType res = Func(f);                                                 \
    reduce(res, rFunc##Op<Type>(), Pstream::msgType(), comm);                 \
    return res;                                                               \
}                                                                             \
TMP_UNARY_FUNCTION(ReturnType, gFunc)

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(Type, gMaxMagSqr, maxMagSqr, maxMagSqr)
G_UNARY_FUNCTION(Type, gMinMagSqr, minMagSqr, minMagSqr)
G_UNARY_FUNCTION(scalar, gSumSqr, sumSqr, sum)
G_UNARY_FUNCTION(scalar, gSumMag, sumMag, sum)
G_UNARY_FUNCTION(Type, gSumCmptMag, sumCmptMag, sum)

#undef G_UNARY_FUNCTION

template<class Type>
scalar gSumProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const int comm
)
{
    scalar SumProd = sumProd(f1, f2);
    reduce(SumProd, sumOp<scalar>(), Pstream::msgType(), comm);
    return SumProd;
}

template<class Type>
Type gSumCmptProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const int comm
)
{
    Type SumProd = sumCmptProd(f1, f2);
    reduce(SumProd, sumOp<Type>(), Pstream::msgType(), comm);
    return SumProd;
}

template<class Type>
Type gAverage
(
    const UList<Type>& f,
    const int comm
)
{
    label n = f.size();
    Type s = sum(f);
    sumReduce(s, n, Pstream::msgType(), comm);

    if (n > 0)
    {
        Type avrg = s/n;

        return avrg;
    }
    else
    {
        WarningIn("gAverage(const UList<Type>&)")
            << "empty field, returning zero." << endl;

        return pTraits<Type>::zero;
    }
}

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

BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
BINARY_OPERATOR(Type, Type, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                 \
                                                                              \
template<class Type1, class Type2>                                            \
void OpFunc                                                                   \
(                                                                             \
    Field<typename product<Type1, Type2>::type>& res,                         \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    TFOR_ALL_F_OP_F_OP_F(productType, res, =, Type1, f1, Op, Type2, f2)       \
}                                                                             \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const UList<Type1>& f1, const UList<Type2>& f2)                   \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<Field<productType> > tRes(new Field<productType>(f1.size()));         \
    OpFunc(tRes(), f1, f2);                                                   \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const UList<Type1>& f1, const tmp<Field<Type2> >& tf2)            \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<Field<productType> > tRes = reuseTmp<productType, Type2>::New(tf2);   \
    OpFunc(tRes(), f1, tf2());                                                \
    reuseTmp<productType, Type2>::clear(tf2);                                 \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const tmp<Field<Type1> >& tf1, const UList<Type2>& f2)            \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<Field<productType> > tRes = reuseTmp<productType, Type1>::New(tf1);   \
    OpFunc(tRes(), tf1(), f2);                                                \
    reuseTmp<productType, Type1>::clear(tf1);                                 \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const tmp<Field<Type1> >& tf1, const tmp<Field<Type2> >& tf2)     \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<Field<productType> > tRes =                                           \
        reuseTmpTmp<productType, Type1, Type1, Type2>::New(tf1, tf2);         \
    OpFunc(tRes(), tf1(), tf2());                                             \
    reuseTmpTmp<productType, Type1, Type1, Type2>::clear(tf1, tf2);           \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
void OpFunc                                                                   \
(                                                                             \
    Field<typename product<Type, Form>::type>& res,                           \
    const UList<Type>& f1,                                                    \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
)                                                                             \
{                                                                             \
    typedef typename product<Type, Form>::type productType;                   \
    TFOR_ALL_F_OP_F_OP_S                                                      \
        (productType, res, =,Type, f1, Op, Form, static_cast<const Form&>(vs))\
}                                                                             \
                                                                              \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
tmp<Field<typename product<Type, Form>::type> >                               \
operator Op(const UList<Type>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs)    \
{                                                                             \
    typedef typename product<Type, Form>::type productType;                   \
    tmp<Field<productType> > tRes(new Field<productType>(f1.size()));         \
    OpFunc(tRes(), f1, static_cast<const Form&>(vs));                         \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
tmp<Field<typename product<Type, Form>::type> >                               \
operator Op                                                                   \
(                                                                             \
    const tmp<Field<Type> >& tf1,                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
)                                                                             \
{                                                                             \
    typedef typename product<Type, Form>::type productType;                   \
    tmp<Field<productType> > tRes = reuseTmp<productType, Type>::New(tf1);    \
    OpFunc(tRes(), tf1(), static_cast<const Form&>(vs));                      \
    reuseTmp<productType, Type>::clear(tf1);                                  \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
void OpFunc                                                                   \
(                                                                             \
    Field<typename product<Form, Type>::type>& res,                           \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const UList<Type>& f1                                                     \
)                                                                             \
{                                                                             \
    typedef typename product<Form, Type>::type productType;                   \
    TFOR_ALL_F_OP_S_OP_F                                                      \
        (productType, res, =,Form,static_cast<const Form&>(vs), Op, Type, f1) \
}                                                                             \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
tmp<Field<typename product<Form, Type>::type> >                               \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const UList<Type>& f1)    \
{                                                                             \
    typedef typename product<Form, Type>::type productType;                   \
    tmp<Field<productType> > tRes(new Field<productType>(f1.size()));         \
    OpFunc(tRes(), static_cast<const Form&>(vs), f1);                         \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
tmp<Field<typename product<Form, Type>::type> >                               \
operator Op                                                                   \
(                                                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs, const tmp<Field<Type> >& tf1      \
)                                                                             \
{                                                                             \
    typedef typename product<Form, Type>::type productType;                   \
    tmp<Field<productType> > tRes = reuseTmp<productType, Type>::New(tf1);    \
    OpFunc(tRes(), static_cast<const Form&>(vs), tf1());                      \
    reuseTmp<productType, Type>::clear(tf1);                                  \
    return tRes;                                                              \
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
