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

#include "GeometricFieldReuseFunctions.H"

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
#define TEMPLATE3                                                              \
    template                                                                   \
    <                                                                          \
        class Type,                                                            \
        class GeoMesh,                                                         \
        template<class> class PrimitiveField1,                                 \
        template<class> class PrimitiveField2,                                 \
        template<class> class PrimitiveField3                                  \
    >
#include "GeometricFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

TEMPLATE2
void component
(
    GeometricField
    <
        typename GeometricField<Type, GeoMesh, PrimitiveField1>::cmptType,
        GeoMesh,
        PrimitiveField1
    >& gcf,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf,
    const direction d
)
{
    component(gcf.primitiveFieldRef(), gf.primitiveField(), d);
    component(gcf.boundaryFieldRef(), gf.boundaryField(), d);
}


TEMPLATE2
void T
(
     GeometricField<Type, GeoMesh, PrimitiveField1>& gf,
     const GeometricField<Type, GeoMesh, PrimitiveField2>& gf1
)
{
    T(gf.primitiveFieldRef(), gf1.primitiveField());
    T(gf.boundaryFieldRef(), gf1.boundaryField());
}


template
<
    class Type,
    class GeoMesh,
    template<class> class PrimitiveField1,
    template<class> class PrimitiveField2,
    direction r
>
void pow
(
    GeometricField
    <
        typename powProduct<Type, r>::type,
        GeoMesh,
        PrimitiveField1
    >& gf,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf1
)
{
    pow(gf.primitiveFieldRef(), gf1.primitiveField(), r);
    pow(gf.boundaryFieldRef(), gf1.boundaryField(), r);
}


template
<
    class Type,
    class GeoMesh,
    template<class> class PrimitiveField,
    direction r
>
tmp<GeometricField<typename powProduct<Type, r>::type, GeoMesh, Field>>
pow
(
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;

    tmp<GeometricField<powProductType, GeoMesh, Field>> tPow
    (
        GeometricField<powProductType, GeoMesh, Field>::New
        (
            "pow(" + gf.name() + ',' + name(r) + ')',
            gf.mesh(),
            pow(gf.dimensions(), r)
        )
    );

    pow<Type, GeoMesh, Field, PrimitiveField, r>(tPow.ref(), gf);

    return tPow;
}


template
<
    class Type,
    class GeoMesh,
    template<class> class PrimitiveField,
    direction r
>
tmp<GeometricField<typename powProduct<Type, r>::type, GeoMesh, Field>>
pow
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;

    const GeometricField<Type, GeoMesh, PrimitiveField>& gf = tgf();

    tmp<GeometricField<powProductType, GeoMesh, Field>> tPow
    (
        GeometricField<powProductType, GeoMesh, Field>::New
        (
            "pow(" + gf.name() + ',' + name(r) + ')',
            gf.mesh(),
            pow(gf.dimensions(), r)
        )
    );

    pow<Type, GeoMesh, Field, PrimitiveField, r>
    (
        tPow.ref(),
        gf
    );

    tgf.clear();

    return tPow;
}


TEMPLATE2
void sqr
(
    GeometricField
    <
        typename outerProduct<Type, Type>::type,
        GeoMesh,
        PrimitiveField1
    >& gf,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf1
)
{
    sqr(gf.primitiveFieldRef(), gf1.primitiveField());
    sqr(gf.boundaryFieldRef(), gf1.boundaryField());
}


TEMPLATE
tmp<GeometricField<typename outerProduct<Type, Type>::type, GeoMesh, Field>>
sqr(const GeometricField<Type, GeoMesh, PrimitiveField>& gf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;

    tmp<GeometricField<outerProductType, GeoMesh, Field>> tSqr
    (
        GeometricField<outerProductType, GeoMesh, Field>::New
        (
            "sqr(" + gf.name() + ')',
            gf.mesh(),
            sqr(gf.dimensions())
        )
    );

    sqr(tSqr.ref(), gf);

    return tSqr;
}


TEMPLATE
tmp<GeometricField<typename outerProduct<Type, Type>::type, GeoMesh, Field>>
sqr(const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;

    const GeometricField<Type, GeoMesh, PrimitiveField>& gf = tgf();

    tmp<GeometricField<outerProductType, GeoMesh, Field>> tSqr
    (
        GeometricField<outerProductType, GeoMesh, Field>::New
        (
            "sqr(" + gf.name() + ')',
            gf.mesh(),
            sqr(gf.dimensions())
        )
    );

    sqr(tSqr.ref(), gf);

    tgf.clear();

    return tSqr;
}


TEMPLATE2
void magSqr
(
    GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf
)
{
    magSqr(gsf.primitiveFieldRef(), gf.primitiveField());
    magSqr(gsf.boundaryFieldRef(), gf.boundaryField());
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> magSqr
(
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf
)
{
    tmp<GeometricField<scalar, GeoMesh, Field>> tMagSqr
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "magSqr(" + gf.name() + ')',
            gf.mesh(),
            sqr(gf.dimensions())
        )
    );

    magSqr(tMagSqr.ref(), gf);

    return tMagSqr;
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> magSqr
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
)
{
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf = tgf();

    tmp<GeometricField<scalar, GeoMesh, Field>> tMagSqr
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "magSqr(" + gf.name() + ')',
            gf.mesh(),
            sqr(gf.dimensions())
        )
    );

    magSqr(tMagSqr.ref(), gf);

    tgf.clear();

    return tMagSqr;
}


TEMPLATE2
void mag
(
    GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf
)
{
    mag(gsf.primitiveFieldRef(), gf.primitiveField());
    mag(gsf.boundaryFieldRef(), gf.boundaryField());
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> mag
(
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf
)
{
    tmp<GeometricField<scalar, GeoMesh, Field>> tMag
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "mag(" + gf.name() + ')',
            gf.mesh(),
            gf.dimensions()
        )
    );

    mag(tMag.ref(), gf);

    return tMag;
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> mag
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
)
{
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf = tgf();

    tmp<GeometricField<scalar, GeoMesh, Field>> tMag
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "mag(" + gf.name() + ')',
            gf.mesh(),
            gf.dimensions()
        )
    );

    mag(tMag.ref(), gf);

    tgf.clear();

    return tMag;
}


TEMPLATE2
void cmptAv
(
    GeometricField
    <
        typename GeometricField<Type, GeoMesh, PrimitiveField1>::cmptType,
        GeoMesh,
        PrimitiveField1
    >& gcf,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf
)
{
    cmptAv(gcf.primitiveFieldRef(), gf.primitiveField());
    cmptAv(gcf.boundaryFieldRef(), gf.boundaryField());
}


TEMPLATE
tmp
<
    GeometricField
    <
        typename GeometricField<Type, GeoMesh, PrimitiveField>::cmptType,
        GeoMesh,
        Field
    >
>
cmptAv(const GeometricField<Type, GeoMesh, PrimitiveField>& gf)
{
    typedef typename
        GeometricField<Type, GeoMesh, PrimitiveField>::cmptType cmptType;

    tmp<GeometricField<cmptType, GeoMesh, Field>> CmptAv
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "cmptAv(" + gf.name() + ')',
            gf.mesh(),
            gf.dimensions()
        )
    );

    cmptAv(CmptAv.ref(), gf);

    return CmptAv;
}


TEMPLATE
tmp
<
    GeometricField
    <
        typename GeometricField<Type, GeoMesh, PrimitiveField>::cmptType,
        GeoMesh,
        Field
    >
>
cmptAv(const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf)
{
    typedef typename
        GeometricField<Type, GeoMesh, PrimitiveField>::cmptType cmptType;

    const GeometricField<Type, GeoMesh, PrimitiveField>& gf = tgf();

    tmp<GeometricField<cmptType, GeoMesh, Field>> CmptAv
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "cmptAv(" + gf.name() + ')',
            gf.mesh(),
            gf.dimensions()
        )
    );

    cmptAv(CmptAv.ref(), gf);

    tgf.clear();

    return CmptAv;
}

UNARY_FUNCTION(Type, Type, cmptMag, cmptMag);


#define UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(returnType, func, gFunc)        \
                                                                               \
TEMPLATE                                                                       \
dimensioned<returnType> func                                                   \
(                                                                              \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf                    \
)                                                                              \
{                                                                              \
    return dimensioned<Type>                                                   \
    (                                                                          \
        #func "(" + gf.name() + ')',                                           \
        gf.dimensions(),                                                       \
        Foam::func(gFunc(gf.primitiveField()), gFunc(gf.boundaryField()))      \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
dimensioned<returnType> func                                                   \
(                                                                              \
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf1             \
)                                                                              \
{                                                                              \
    dimensioned<returnType> res = func(tgf1());                                \
    tgf1.clear();                                                              \
    return res;                                                                \
}

UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(Type, max, gMax)
UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(Type, min, gMin)

#undef UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY


#define UNARY_REDUCTION_FUNCTION(returnType, func, gFunc)                      \
                                                                               \
TEMPLATE                                                                       \
dimensioned<returnType> func                                                   \
(                                                                              \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf                    \
)                                                                              \
{                                                                              \
    return dimensioned<Type>                                                   \
    (                                                                          \
        #func "(" + gf.name() + ')',                                           \
        gf.dimensions(),                                                       \
        gFunc(gf.primitiveField())                                             \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
dimensioned<returnType> func                                                   \
(                                                                              \
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf1             \
)                                                                              \
{                                                                              \
    dimensioned<returnType> res = func(tgf1());                                \
    tgf1.clear();                                                              \
    return res;                                                                \
}

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
    template<class> class PrimitiveField2,                                     \
    template<class> class PrimitiveField3                                      \
>                                                                              \
void opFunc                                                                    \
(                                                                              \
    GeometricField                                                             \
    <                                                                          \
        typename product<Type1, Type2>::type,                                  \
        GeoMesh,                                                               \
        PrimitiveField1                                                        \
    >& gf,                                                                     \
    const GeometricField<Type1, GeoMesh, PrimitiveField2>& gf1,                \
    const GeometricField<Type2, GeoMesh, PrimitiveField3>& gf2                 \
)                                                                              \
{                                                                              \
    Foam::opFunc                                                               \
    (                                                                          \
        gf.primitiveFieldRef(),                                                \
        gf1.primitiveField(),                                                  \
        gf2.primitiveField()                                                   \
    );                                                                         \
    Foam::opFunc                                                               \
    (                                                                          \
        gf.boundaryFieldRef(),                                                 \
        gf1.boundaryField(),                                                   \
        gf2.boundaryField()                                                    \
    );                                                                         \
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
tmp<GeometricField<typename product<Type1, Type2>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField1>& gf1,                \
    const GeometricField<Type2, GeoMesh, PrimitiveField2>& gf2                 \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    tmp<GeometricField<productType, GeoMesh, Field>> tRes                      \
    (                                                                          \
        GeometricField<productType, GeoMesh, Field>::New                       \
        (                                                                      \
            '(' + gf1.name() + #op + gf2.name() + ')',                         \
            gf1.mesh(),                                                        \
            gf1.dimensions() op gf2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::opFunc(tRes.ref(), gf1, gf2);                                        \
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
tmp<GeometricField<typename product<Type1, Type2>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField1>& gf1,                \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField2>>& tgf2           \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    const GeometricField<Type2, GeoMesh, PrimitiveField2>& gf2 = tgf2();       \
                                                                               \
    tmp<GeometricField<productType, GeoMesh, Field>> tRes                      \
    (                                                                          \
        reuseTmpGeometricField                                                 \
        <productType, Type2, GeoMesh, PrimitiveField2>::New                    \
        (                                                                      \
            tgf2,                                                              \
            '(' + gf1.name() + #op + gf2.name() + ')',                         \
            gf1.dimensions() op gf2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::opFunc(tRes.ref(), gf1, gf2);                                        \
                                                                               \
    tgf2.clear();                                                              \
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
tmp<GeometricField<typename product<Type1, Type2>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField1>>& tgf1,          \
    const GeometricField<Type2, GeoMesh, PrimitiveField2>& gf2                 \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    const GeometricField<Type1, GeoMesh, PrimitiveField1>& gf1 = tgf1();       \
                                                                               \
    tmp<GeometricField<productType, GeoMesh, Field>> tRes                      \
    (                                                                          \
        reuseTmpGeometricField                                                 \
        <productType, Type1, GeoMesh, PrimitiveField1>::New                    \
        (                                                                      \
            tgf1,                                                              \
            '(' + gf1.name() + #op + gf2.name() + ')',                         \
            gf1.dimensions() op gf2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::opFunc(tRes.ref(), gf1, gf2);                                        \
                                                                               \
    tgf1.clear();                                                              \
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
tmp<GeometricField<typename product<Type1, Type2>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField1>>& tgf1,          \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField2>>& tgf2           \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    const GeometricField<Type1, GeoMesh, PrimitiveField1>& gf1 = tgf1();       \
    const GeometricField<Type2, GeoMesh, PrimitiveField2>& gf2 = tgf2();       \
                                                                               \
    tmp<GeometricField<productType, GeoMesh, Field>> tRes                      \
    (                                                                          \
        reuseTmpTmpGeometricField                                              \
        <                                                                      \
            productType,                                                       \
            Type1,                                                             \
            Type2,                                                             \
            GeoMesh,                                                           \
            PrimitiveField1,                                                   \
            PrimitiveField2                                                    \
        >::New                                                                 \
        (                                                                      \
            tgf1,                                                              \
            tgf2,                                                              \
            '(' + gf1.name() + #op + gf2.name() + ')',                         \
            gf1.dimensions() op gf2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::opFunc(tRes.ref(), gf1, gf2);                                        \
                                                                               \
    tgf1.clear();                                                              \
    tgf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField1,                                     \
    template<class> class PrimitiveField2                                      \
>                                                                              \
void opFunc                                                                    \
(                                                                              \
    GeometricField                                                             \
    <                                                                          \
        typename product<Type, Form>::type,                                    \
        GeoMesh,                                                               \
        PrimitiveField1                                                        \
    >& gf,                                                                     \
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf1,                 \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    Foam::opFunc(gf.primitiveFieldRef(), gf1.primitiveField(), dvs.value());   \
    Foam::opFunc(gf.boundaryFieldRef(), gf1.boundaryField(), dvs.value());     \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<GeometricField<typename product<Form, Type>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf1,                  \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
                                                                               \
    tmp<GeometricField<productType, GeoMesh, Field>> tRes                      \
    (                                                                          \
        GeometricField<productType, GeoMesh, Field>::New                       \
        (                                                                      \
            '(' + gf1.name() + #op + dvs.name() + ')',                         \
            gf1.mesh(),                                                        \
            gf1.dimensions() op dvs.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::opFunc(tRes.ref(), gf1, dvs);                                        \
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
tmp<GeometricField<typename product<Form, Type>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf1,                  \
    const VectorSpace<Form, Cmpt, nCmpt>& vs                                   \
)                                                                              \
{                                                                              \
    return gf1 op dimensioned<Form>(static_cast<const Form&>(vs));             \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<GeometricField<typename product<Type, Form>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf1,            \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
                                                                               \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf1 = tgf1();         \
                                                                               \
    tmp<GeometricField<productType, GeoMesh, Field>> tRes                      \
    (                                                                          \
        reuseTmpGeometricField                                                 \
        <productType, Type, GeoMesh, PrimitiveField>::New                      \
        (                                                                      \
            tgf1,                                                              \
            '(' + gf1.name() + #op + dvs.name() + ')',                         \
            gf1.dimensions() op dvs.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::opFunc(tRes.ref(), gf1, dvs);                                        \
                                                                               \
    tgf1.clear();                                                              \
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
tmp<GeometricField<typename product<Form, Type>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf1,            \
    const VectorSpace<Form, Cmpt, nCmpt>& vs                                   \
)                                                                              \
{                                                                              \
    return tgf1 op dimensioned<Form>(static_cast<const Form&>(vs));            \
}                                                                              \
                                                                               \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField1,                                     \
    template<class> class PrimitiveField2                                      \
>                                                                              \
void opFunc                                                                    \
(                                                                              \
    GeometricField                                                             \
    <                                                                          \
        typename product<Form, Type>::type,                                    \
        GeoMesh,                                                               \
        PrimitiveField1                                                        \
    >& gf,                                                                     \
    const dimensioned<Form>& dvs,                                              \
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf1                  \
)                                                                              \
{                                                                              \
    Foam::opFunc(gf.primitiveFieldRef(), dvs.value(), gf1.primitiveField());   \
    Foam::opFunc(gf.boundaryFieldRef(), dvs.value(), gf1.boundaryField());     \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<GeometricField<typename product<Form, Type>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf1                   \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
                                                                               \
    tmp<GeometricField<productType, GeoMesh, Field>> tRes                      \
    (                                                                          \
        GeometricField<productType, GeoMesh, Field>::New                       \
        (                                                                      \
            '(' + dvs.name() + #op + gf1.name() + ')',                         \
            gf1.mesh(),                                                        \
            dvs.dimensions() op gf1.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::opFunc(tRes.ref(), dvs, gf1);                                        \
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
tmp<GeometricField<typename product<Form, Type>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const VectorSpace<Form, Cmpt, nCmpt>& vs,                                  \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf1                   \
)                                                                              \
{                                                                              \
    return dimensioned<Form>(static_cast<const Form&>(vs)) op gf1;             \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Type,                                                                \
    class GeoMesh,                                                             \
    template<class> class PrimitiveField                                       \
>                                                                              \
tmp<GeometricField<typename product<Form, Type>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf1             \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
                                                                               \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf1 = tgf1();         \
                                                                               \
    tmp<GeometricField<productType, GeoMesh, Field>> tRes                      \
    (                                                                          \
        reuseTmpGeometricField                                                 \
        <productType, Type, GeoMesh, PrimitiveField>::New                      \
        (                                                                      \
            tgf1,                                                              \
            '(' + dvs.name() + #op + gf1.name() + ')',                         \
            dvs.dimensions() op gf1.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::opFunc(tRes.ref(), dvs, gf1);                                        \
                                                                               \
    tgf1.clear();                                                              \
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
tmp<GeometricField<typename product<Form, Type>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const VectorSpace<Form, Cmpt, nCmpt>& vs,                                  \
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf1             \
)                                                                              \
{                                                                              \
    return dimensioned<Form>(static_cast<const Form&>(vs)) op tgf1;            \
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
