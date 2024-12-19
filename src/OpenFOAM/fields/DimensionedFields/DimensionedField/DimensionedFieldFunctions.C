/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "DimensionedFieldReuseFunctions.H"

#define TEMPLATE                                                               \
    template                                                                   \
    <                                                                          \
        class Type,                                                            \
        class GeoMesh,                                                         \
        template<class> class PrimitiveField                                   \
    >
#define TEMPLATE2                                                              \
    template                                                                   \
    <                                                                          \
        class Type,                                                            \
        class GeoMesh,                                                         \
        template<class> class PrimitiveField1,                                 \
        template<class> class PrimitiveField2                                  \
    >
#include "DimensionedFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template
<
    class Type,
    class GeoMesh,
    template<class> class PrimitiveField,
    direction r
>
tmp<DimensionedField<typename powProduct<Type, r>::type, GeoMesh, Field>>
pow
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;

    tmp<DimensionedField<powProductType, GeoMesh, Field>> tPow
    (
        DimensionedField<powProductType, GeoMesh, Field>::New
        (
            "pow(" + df.name() + ',' + name(r) + ')',
            df.mesh(),
            pow(df.dimensions(), r)
        )
    );

    pow<Type, r, GeoMesh>(tPow.ref().primitiveFieldRef(), df.primitiveField());

    return tPow;
}


template
<
    class Type,
    class GeoMesh,
    template<class> class PrimitiveField,
    direction r
>
tmp<DimensionedField<typename powProduct<Type, r>::type, GeoMesh, Field>>
pow
(
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;

    const DimensionedField<Type, GeoMesh, PrimitiveField>& df = tdf();

    tmp<DimensionedField<powProductType, GeoMesh, Field>> tPow =
        reuseTmpDimensionedField
        <
            powProductType,
            Type,
            GeoMesh,
            PrimitiveField
        >::New
        (
            tdf,
            "pow(" + df.name() + ',' + name(r) + ')',
            pow(df.dimensions(), r)
        );

    pow<Type, r, GeoMesh>(tPow.ref().primitiveFieldRef(), df.primitiveField());

    tdf.clear();

    return tPow;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp
<
    DimensionedField
    <
        typename outerProduct<Type, Type>::type,
        GeoMesh,
        Field
    >
>
sqr(const DimensionedField<Type, GeoMesh, PrimitiveField>& df)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;

    tmp<DimensionedField<outerProductType, GeoMesh, Field>> tSqr
    (
        DimensionedField<outerProductType, GeoMesh, Field>::New
        (
            "sqr(" + df.name() + ')',
            df.mesh(),
            sqr(df.dimensions())
        )
    );

    sqr(tSqr.ref().primitiveFieldRef(), df.primitiveField());

    return tSqr;
}

template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp
<
    DimensionedField
    <
        typename outerProduct<Type, Type>::type,
        GeoMesh,
        Field
    >
>
sqr(const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;

    const DimensionedField<Type, GeoMesh, PrimitiveField>& df = tdf();

    tmp<DimensionedField<outerProductType, GeoMesh, Field>> tSqr =
        reuseTmpDimensionedField
        <
            outerProductType,
            Type,
            GeoMesh,
            PrimitiveField
        >::New
        (
            tdf,
            "sqr(" + df.name() + ')',
            sqr(df.dimensions())
        );

    sqr(tSqr.ref().primitiveFieldRef(), df.primitiveField());

    tdf.clear();

    return tSqr;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> magSqr
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df
)
{
    tmp<DimensionedField<scalar, GeoMesh, Field>> tMagSqr
    (
        DimensionedField<scalar, GeoMesh, Field>::New
        (
            "magSqr(" + df.name() + ')',
            df.mesh(),
            sqr(df.dimensions())
        )
    );

    magSqr(tMagSqr.ref().primitiveFieldRef(), df.primitiveField());

    return tMagSqr;
}

template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> magSqr
(
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf
)
{
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df = tdf();

    tmp<DimensionedField<scalar, GeoMesh, Field>> tMagSqr =
        reuseTmpDimensionedField<scalar, Type, GeoMesh, PrimitiveField>::New
        (
            tdf,
            "magSqr(" + df.name() + ')',
            sqr(df.dimensions())
        );

    magSqr(tMagSqr.ref().primitiveFieldRef(), df.primitiveField());

    tdf.clear();

    return tMagSqr;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> mag
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df
)
{
    tmp<DimensionedField<scalar, GeoMesh, Field>> tMag
    (
        DimensionedField<scalar, GeoMesh, Field>::New
        (
            "mag(" + df.name() + ')',
            df.mesh(),
            df.dimensions()
        )
    );

    mag(tMag.ref().primitiveFieldRef(), df.primitiveField());

    return tMag;
}

template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> mag
(
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf
)
{
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df = tdf();

    tmp<DimensionedField<scalar, GeoMesh, Field>> tMag =
        reuseTmpDimensionedField<scalar, Type, GeoMesh, PrimitiveField>::New
        (
            tdf,
            "mag(" + df.name() + ')',
            df.dimensions()
        );

    mag(tMag.ref().primitiveFieldRef(), df.primitiveField());

    tdf.clear();

    return tMag;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp
<
    DimensionedField
    <
        typename DimensionedField<Type, GeoMesh, PrimitiveField>::cmptType,
        GeoMesh,
        Field
    >
>
cmptAv(const DimensionedField<Type, GeoMesh, PrimitiveField>& df)
{
    typedef typename
        DimensionedField<Type, GeoMesh, PrimitiveField>::cmptType cmptType;

    tmp<DimensionedField<cmptType, GeoMesh, Field>> CmptAv
    (
        DimensionedField<scalar, GeoMesh, PrimitiveField>::New
        (
            "cmptAv(" + df.name() + ')',
            df.mesh(),
            df.dimensions()
        )
    );

    cmptAv(CmptAv.ref().primitiveFieldRef(), df.primitiveField());

    return CmptAv;
}

template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp
<
    DimensionedField
    <
        typename DimensionedField<Type, GeoMesh, PrimitiveField>::cmptType,
        GeoMesh,
        Field
    >
>
cmptAv(const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf)
{
    typedef typename DimensionedField<Type, GeoMesh, Field>::cmptType
        cmptType;

    const DimensionedField<Type, GeoMesh, PrimitiveField>& df = tdf();

    tmp<DimensionedField<cmptType, GeoMesh, Field>> CmptAv =
        reuseTmpDimensionedField<cmptType, Type, GeoMesh, PrimitiveField>::New
        (
            tdf,
            "cmptAv(" + df.name() + ')',
            df.dimensions()
        );

    cmptAv(CmptAv.ref().primitiveFieldRef(), df.primitiveField());

    tdf.clear();

    return CmptAv;
}

UNARY_FUNCTION(Type, Type, cmptMag, cmptMag);


#define UNARY_REDUCTION_FUNCTION(returnType, func, dfunc)                      \
                                                                               \
template<class Type, class GeoMesh, template<class> class PrimitiveField>      \
dimensioned<returnType> func                                                   \
(                                                                              \
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df                  \
)                                                                              \
{                                                                              \
    return dimensioned<Type>                                                   \
    (                                                                          \
        #func "(" + df.name() + ')',                                           \
        df.dimensions(),                                                       \
        dfunc(df.primitiveField())                                             \
    );                                                                         \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh, template<class> class PrimitiveField>      \
dimensioned<returnType> func                                                   \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf1           \
)                                                                              \
{                                                                              \
    dimensioned<returnType> res = func(tdf1());                                \
    tdf1.clear();                                                              \
    return res;                                                                \
}

UNARY_REDUCTION_FUNCTION(Type, max, gMax)
UNARY_REDUCTION_FUNCTION(Type, min, gMin)
UNARY_REDUCTION_FUNCTION(Type, sum, gSum)
UNARY_REDUCTION_FUNCTION(scalar, sumMag, gSumMag)
UNARY_REDUCTION_FUNCTION(Type, average, gAverage)

#undef UNARY_REDUCTION_FUNCTION


BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(Type, Type, -, negate, transform)

BINARY_OPERATOR(Type, Type, scalar, *, '*', multiply)
BINARY_OPERATOR(Type, scalar, Type, *, '*', multiply)
BINARY_OPERATOR(Type, Type, scalar, /, '|', divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, '*', multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, '*', multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, '|', divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, op, opFunc)                                  \
                                                                               \
template                                                                       \
<                                                                              \
    class Type1,                                                               \
    class Type2,                                                               \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField1,                                     \
    template<class> class PrimitiveField2                                      \
>                                                                              \
tmp                                                                            \
<                                                                              \
    DimensionedField                                                           \
    <                                                                          \
        typename product<Type1, Type2>::type,                                  \
        GeoMesh,                                                               \
        Field                                                                  \
    >                                                                          \
>                                                                              \
operator op                                                                    \
(                                                                              \
    const DimensionedField<Type1, GeoMesh, PrimitiveField1>& df1,              \
    const DimensionedField<Type2, GeoMesh, PrimitiveField2>& df2               \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<DimensionedField<productType, GeoMesh, Field>> tRes                    \
    (                                                                          \
        DimensionedField<productType, GeoMesh, Field>::New                     \
        (                                                                      \
            '(' + df1.name() + #op + df2.name() + ')',                         \
            df1.mesh(),                                                        \
            df1.dimensions() op df2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::opFunc                                                               \
    (                                                                          \
        tRes.ref().primitiveFieldRef(),                                        \
        df1.primitiveField(),                                                  \
        df2.primitiveField()                                                   \
    );                                                                         \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type1,                                                               \
    class Type2,                                                               \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField1,                                     \
    template<class> class PrimitiveField2                                      \
>                                                                              \
tmp                                                                            \
<                                                                              \
    DimensionedField                                                           \
    <                                                                          \
        typename product<Type1, Type2>::type,                                  \
        GeoMesh,                                                               \
        Field                                                                  \
    >                                                                          \
>                                                                              \
operator op                                                                    \
(                                                                              \
    const DimensionedField<Type1, GeoMesh, PrimitiveField1>& df1,              \
    const tmp<DimensionedField<Type2, GeoMesh, PrimitiveField2>>& tdf2         \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    const DimensionedField<Type2, GeoMesh, PrimitiveField2>& df2 = tdf2();     \
                                                                               \
    tmp<DimensionedField<productType, GeoMesh, Field>> tRes =                  \
        reuseTmpDimensionedField                                               \
        <                                                                      \
            productType,                                                       \
            Type2,                                                             \
            GeoMesh,                                                           \
            PrimitiveField2                                                    \
        >::New                                                                 \
        (                                                                      \
            tdf2,                                                              \
            '(' + df1.name() + #op + df2.name() + ')',                         \
            df1.dimensions() op df2.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc                                                               \
    (                                                                          \
        tRes.ref().primitiveFieldRef(),                                        \
        df1.primitiveField(),                                                  \
        df2.primitiveField()                                                   \
    );                                                                         \
                                                                               \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type1,                                                               \
    class Type2,                                                               \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField1,                                     \
    template<class> class PrimitiveField2                                      \
>                                                                              \
tmp                                                                            \
<                                                                              \
    DimensionedField                                                           \
    <                                                                          \
        typename product<Type1, Type2>::type,                                  \
        GeoMesh,                                                               \
        Field                                                                  \
    >                                                                          \
>                                                                              \
operator op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh, PrimitiveField1>>& tdf1,        \
    const DimensionedField<Type2, GeoMesh, PrimitiveField2>& df2               \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    const DimensionedField<Type1, GeoMesh, PrimitiveField1>& df1 = tdf1();     \
                                                                               \
    tmp<DimensionedField<productType, GeoMesh, Field>> tRes =                  \
        reuseTmpDimensionedField                                               \
        <                                                                      \
            productType,                                                       \
            Type1,                                                             \
            GeoMesh,                                                           \
            PrimitiveField1                                                    \
        >::New                                                                 \
        (                                                                      \
            tdf1,                                                              \
            '(' + df1.name() + #op + df2.name() + ')',                         \
            df1.dimensions() op df2.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc                                                               \
    (                                                                          \
        tRes.ref().primitiveFieldRef(),                                        \
        df1.primitiveField(),                                                  \
        df2.primitiveField()                                                   \
    );                                                                         \
                                                                               \
    tdf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type1,                                                               \
    class Type2,                                                               \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField1,                                     \
    template<class> class PrimitiveField2                                      \
>                                                                              \
tmp                                                                            \
<                                                                              \
    DimensionedField                                                           \
    <                                                                          \
        typename product<Type1, Type2>::type,                                  \
        GeoMesh,                                                               \
        Field                                                                  \
    >                                                                          \
>                                                                              \
operator op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh, PrimitiveField1>>& tdf1,        \
    const tmp<DimensionedField<Type2, GeoMesh, PrimitiveField2>>& tdf2         \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    const DimensionedField<Type1, GeoMesh, PrimitiveField1>& df1 = tdf1();     \
    const DimensionedField<Type2, GeoMesh, PrimitiveField2>& df2 = tdf2();     \
                                                                               \
    tmp<DimensionedField<productType, GeoMesh>> tRes =                         \
        reuseTmpTmpDimensionedField                                            \
        <                                                                      \
            productType,                                                       \
            Type1,                                                             \
            Type2,                                                             \
            GeoMesh,                                                           \
            PrimitiveField1,                                                   \
            PrimitiveField2                                                    \
        >::New                                                                 \
        (                                                                      \
            tdf1,                                                              \
            tdf2,                                                              \
            '(' + df1.name() + #op + df2.name() + ')',                         \
            df1.dimensions() op df2.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc                                                               \
    (                                                                          \
        tRes.ref().primitiveFieldRef(),                                        \
        df1.primitiveField(),                                                  \
        df2.primitiveField()                                                   \
    );                                                                         \
                                                                               \
    tdf1.clear();                                                              \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df1,                \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
                                                                               \
    tmp<DimensionedField<productType, GeoMesh, Field>> tRes                    \
    (                                                                          \
        DimensionedField<productType, GeoMesh, Field>::New                     \
        (                                                                      \
            '(' + df1.name() + #op + dvs.name() + ')',                         \
            df1.mesh(),                                                        \
            df1.dimensions() op dvs.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::opFunc                                                               \
    (                                                                          \
        tRes.ref().primitiveFieldRef(),                                        \
        df1.primitiveField(),                                                  \
        dvs.value()                                                            \
    );                                                                         \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df1,                \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    return df1 op dimensioned<Form>(static_cast<const Form&>(vs));             \
}                                                                              \
                                                                               \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf1,          \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
                                                                               \
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df1 = tdf1();       \
                                                                               \
    tmp<DimensionedField<productType, GeoMesh, Field>> tRes =                  \
        reuseTmpDimensionedField                                               \
        <                                                                      \
            productType,                                                       \
            Type,                                                              \
            GeoMesh,                                                           \
            PrimitiveField                                                     \
        >::New                                                                 \
        (                                                                      \
            tdf1,                                                              \
            '(' + df1.name() + #op + dvs.name() + ')',                         \
            df1.dimensions() op dvs.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc                                                               \
    (                                                                          \
        tRes.ref().primitiveFieldRef(),                                        \
        df1.primitiveField(),                                                  \
        dvs.value()                                                            \
    );                                                                         \
                                                                               \
    tdf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf1,          \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    return tdf1 op dimensioned<Form>(static_cast<const Form&>(vs));            \
}                                                                              \
                                                                               \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df1                 \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<DimensionedField<productType, GeoMesh, Field>> tRes                    \
    (                                                                          \
        DimensionedField<productType, GeoMesh, Field>::New                     \
        (                                                                      \
            '(' + dvs.name() + #op + df1.name() + ')',                         \
            df1.mesh(),                                                        \
            dvs.dimensions() op df1.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::opFunc                                                               \
    (                                                                          \
        tRes.ref().primitiveFieldRef(),                                        \
        dvs.value(),                                                           \
        df1.primitiveField()                                                   \
    );                                                                         \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df1                 \
)                                                                              \
{                                                                              \
    return dimensioned<Form>(static_cast<const Form&>(vs)) op df1;             \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf1           \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
                                                                               \
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df1 = tdf1();       \
                                                                               \
    tmp<DimensionedField<productType, GeoMesh, Field>> tRes =                  \
        reuseTmpDimensionedField                                               \
        <                                                                      \
            productType,                                                       \
            Type,                                                              \
            GeoMesh,                                                           \
            PrimitiveField                                                     \
        >::New                                                                 \
        (                                                                      \
            tdf1,                                                              \
            '(' + dvs.name() + #op + df1.name() + ')',                         \
            dvs.dimensions() op df1.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc                                                               \
    (                                                                          \
        tRes.ref().primitiveFieldRef(),                                        \
        dvs.value(),                                                           \
        df1.primitiveField()                                                   \
    );                                                                         \
                                                                               \
    tdf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf1           \
)                                                                              \
{                                                                              \
    return dimensioned<Form>(static_cast<const Form&>(vs)) op tdf1;            \
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
