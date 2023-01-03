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

#include "DimensionedScalarField.H"

#define TEMPLATE template<class GeoMesh, template<class> class PrimitiveField>
#define TEMPLATE2                                                              \
    template                                                                   \
    <                                                                          \
        class GeoMesh,                                                         \
        template<class> class PrimitiveField1,                                 \
        template<class> class PrimitiveField2                                  \
    >
#include "DimensionedFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> stabilise
(
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf,
    const dimensioned<scalar>& ds
)
{
    tmp<DimensionedField<scalar, GeoMesh, Field>> tRes
    (
        DimensionedField<scalar, GeoMesh, Field>::New
        (
            "stabilise(" + dsf.name() + ',' + ds.name() + ')',
            dsf.mesh(),
            dsf.dimensions() + ds.dimensions()
        )
    );

    stabilise(tRes.ref().primitiveFieldRef(), dsf.primitiveField(), ds.value());

    return tRes;
}


template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> stabilise
(
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField>>& tdsf,
    const dimensioned<scalar>& ds
)
{
    const DimensionedField<scalar, GeoMesh, Field>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh, Field>> tRes = New
    (
        tdsf,
        "stabilise(" + dsf.name() + ',' + ds.name() + ')',
        dsf.dimensions() + ds.dimensions()
    );

    stabilise(tRes.ref().primitiveFieldRef(), dsf.primitiveField(), ds.value());

    tdsf.clear();

    return tRes;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(scalar, scalar, scalar, +, '+', add)
BINARY_TYPE_OPERATOR(scalar, scalar, scalar, -, '-', subtract)

BINARY_OPERATOR(scalar, scalar, scalar, *, '*', multiply)
BINARY_OPERATOR(scalar, scalar, scalar, /, '|', divide)

BINARY_TYPE_OPERATOR_SF(scalar, scalar, scalar, /, '|', divide)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf1,
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf2
)
{
    if (!dsf1.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base field is not dimensionless: " << dsf1.dimensions()
            << exit(FatalError);
    }

    if (!dsf2.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << dsf2.dimensions()
            << exit(FatalError);
    }

    tmp<DimensionedField<scalar, GeoMesh, Field>> tPow
    (
        DimensionedField<scalar, GeoMesh, Field>::New
        (
            "pow(" + dsf1.name() + ',' + dsf2.name() + ')',
            dsf1.mesh(),
            dimless
        )
    );

    pow
    (
        tPow.ref().primitiveFieldRef(),
        dsf1.primitiveField(),
        dsf2.primitiveField()
    );

    return tPow;
}


template
<
    class GeoMesh,
    template<class> class PrimitiveField1,
    template<class> class PrimitiveField2
>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField1>>& tdsf1,
    const DimensionedField<scalar, GeoMesh, PrimitiveField2>& dsf2
)
{
    const DimensionedField<scalar, GeoMesh, PrimitiveField1>& dsf1 = tdsf1();

    if (!dsf1.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base field is not dimensionless: " << dsf1.dimensions()
            << exit(FatalError);
    }

    if (!dsf2.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << dsf2.dimensions()
            << exit(FatalError);
    }

    tmp<DimensionedField<scalar, GeoMesh, Field>> tPow = New
    (
        tdsf1,
        "pow(" + dsf1.name() + ',' + dsf2.name() + ')',
        dimless
    );

    pow
    (
        tPow.ref().primitiveFieldRef(),
        dsf1.primitiveField(),
        dsf2.primitiveField()
    );

    tdsf1.clear();

    return tPow;
}


template
<
    class GeoMesh,
    template<class> class PrimitiveField1,
    template<class> class PrimitiveField2
>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const DimensionedField<scalar, GeoMesh, PrimitiveField1>& dsf1,
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField2>>& tdsf2
)
{
    const DimensionedField<scalar, GeoMesh, PrimitiveField2>& dsf2 = tdsf2();

    if (!dsf1.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base field is not dimensionless: " << dsf1.dimensions()
            << exit(FatalError);
    }

    if (!dsf2.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << dsf2.dimensions()
            << exit(FatalError);
    }

    tmp<DimensionedField<scalar, GeoMesh, Field>> tPow = New
    (
        tdsf2,
        "pow(" + dsf1.name() + ',' + dsf2.name() + ')',
        dimless
    );

    pow
    (
        tPow.ref().primitiveFieldRef(),
        dsf1.primitiveField(),
        dsf2.primitiveField()
    );

    tdsf2.clear();

    return tPow;
}


template
<
    class GeoMesh,
    template<class> class PrimitiveField1,
    template<class> class PrimitiveField2
>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField1>>& tdsf1,
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField2>>& tdsf2
)
{
    const DimensionedField<scalar, GeoMesh, PrimitiveField1>& dsf1 = tdsf1();
    const DimensionedField<scalar, GeoMesh, PrimitiveField2>& dsf2 = tdsf2();

    if (!dsf1.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base field is not dimensionless: " << dsf1.dimensions()
            << exit(FatalError);
    }

    if (!dsf2.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << dsf2.dimensions()
            << exit(FatalError);
    }

    tmp<DimensionedField<scalar, GeoMesh>> tPow =
        reuseTmpTmpDimensionedField
        <
            scalar,
            scalar,
            scalar,
            GeoMesh,
            PrimitiveField1,
            PrimitiveField2
        >::New
        (
            tdsf1,
            tdsf2,
            "pow(" + dsf1.name() + ',' + dsf2.name() + ')',
            dimless
        );

    pow
    (
        tPow.ref().primitiveFieldRef(),
        dsf1.primitiveField(),
        dsf2.primitiveField()
    );

    tdsf1.clear();
    tdsf2.clear();

    return tPow;
}


template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf,
    const dimensionedScalar& ds
)
{
    if (!ds.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent is not dimensionless: " << ds.dimensions()
            << exit(FatalError);
    }

    tmp<DimensionedField<scalar, GeoMesh, Field>> tPow
    (
        DimensionedField<scalar, GeoMesh, Field>::New
        (
            "pow(" + dsf.name() + ',' + ds.name() + ')',
            dsf.mesh(),
            pow(dsf.dimensions(), ds)
        )
    );

    pow(tPow.ref().primitiveFieldRef(), dsf.primitiveField(), ds.value());

    return tPow;
}


template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField>>& tdsf,
    const dimensionedScalar& ds
)
{
    if (!ds.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent is not dimensionless: " << ds.dimensions()
            << exit(FatalError);
    }

    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh, Field>> tPow = New
    (
        tdsf,
        "pow(" + dsf.name() + ',' + ds.name() + ')',
        pow(dsf.dimensions(), ds)
    );

    pow(tPow.ref().primitiveFieldRef(), dsf.primitiveField(), ds.value());

    tdsf.clear();

    return tPow;
}


template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf,
    const scalar& s
)
{
    return pow(dsf, dimensionedScalar(s));
}


template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField>>& tdsf,
    const scalar& s
)
{
    return pow(tdsf, dimensionedScalar(s));
}


template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const dimensionedScalar& ds,
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf
)
{
    if (!ds.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base scalar is not dimensionless: " << ds.dimensions()
            << exit(FatalError);
    }

    if (!dsf.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << dsf.dimensions()
            << exit(FatalError);
    }

    tmp<DimensionedField<scalar, GeoMesh, Field>> tPow
    (
        DimensionedField<scalar, GeoMesh, Field>::New
        (
            "pow(" + ds.name() + ',' + dsf.name() + ')',
            dsf.mesh(),
            dimless
        )
    );

    pow(tPow.ref().primitiveFieldRef(), ds.value(), dsf.primitiveField());

    return tPow;
}


template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const dimensionedScalar& ds,
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField>>& tdsf
)
{
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf = tdsf();

    if (!ds.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base scalar is not dimensionless: " << ds.dimensions()
            << exit(FatalError);
    }

    if (!dsf.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << dsf.dimensions()
            << exit(FatalError);
    }

    tmp<DimensionedField<scalar, GeoMesh, Field>> tPow = New
    (
        tdsf,
        "pow(" + ds.name() + ',' + dsf.name() + ')',
        dimless
    );

    pow(tPow.ref().primitiveFieldRef(), ds.value(), dsf.primitiveField());

    tdsf.clear();

    return tPow;
}

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const scalar& s,
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf
)
{
    return pow(dimensionedScalar(s), dsf);
}

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> pow
(
    const scalar& s,
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField>>& tdsf
)
{
    return pow(dimensionedScalar(s), tdsf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template
<
    class GeoMesh,
    template<class> class PrimitiveField1,
    template<class> class PrimitiveField2
>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const DimensionedField<scalar, GeoMesh, PrimitiveField1>& dsf1,
    const DimensionedField<scalar, GeoMesh, PrimitiveField2>& dsf2
)
{
    tmp<DimensionedField<scalar, GeoMesh, Field>> tAtan2
    (
        DimensionedField<scalar, GeoMesh, Field>::New
        (
            "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
            dsf1.mesh(),
            atan2(dsf1.dimensions(), dsf2.dimensions())
        )
    );

    atan2
    (
        tAtan2.ref().primitiveFieldRef(),
        dsf1.primitiveField(),
        dsf2.primitiveField()
    );

    return tAtan2;
}


template
<
    class GeoMesh,
    template<class> class PrimitiveField1,
    template<class> class PrimitiveField2
>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField1>>& tdsf1,
    const DimensionedField<scalar, GeoMesh, PrimitiveField2>& dsf2
)
{
    const DimensionedField<scalar, GeoMesh, PrimitiveField1>& dsf1 = tdsf1();

    tmp<DimensionedField<scalar, GeoMesh, Field>> tAtan2 = New
    (
        tdsf1,
        "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
        atan2(dsf1.dimensions(), dsf2.dimensions())
    );

    atan2
    (
        tAtan2.ref().primitiveFieldRef(),
        dsf1.primitiveField(),
        dsf2.primitiveField()
    );

    tdsf1.clear();

    return tAtan2;
}


template
<
    class GeoMesh,
    template<class> class PrimitiveField1,
    template<class> class PrimitiveField2
>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const DimensionedField<scalar, GeoMesh, PrimitiveField1>& dsf1,
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField2>>& tdsf2
)
{
    const DimensionedField<scalar, GeoMesh, PrimitiveField2>& dsf2 = tdsf2();

    tmp<DimensionedField<scalar, GeoMesh, Field>> tAtan2 = New
    (
        tdsf2,
        "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
        atan2(dsf1.dimensions(), dsf2.dimensions())
    );

    atan2
    (
        tAtan2.ref().primitiveFieldRef(),
        dsf1.primitiveField(),
        dsf2.primitiveField()
    );

    tdsf2.clear();

    return tAtan2;
}

template
<
    class GeoMesh,
    template<class> class PrimitiveField1,
    template<class> class PrimitiveField2
>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField1>>& tdsf1,
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField2>>& tdsf2
)
{
    const DimensionedField<scalar, GeoMesh, PrimitiveField1>& dsf1 = tdsf1();
    const DimensionedField<scalar, GeoMesh, PrimitiveField2>& dsf2 = tdsf2();

    tmp<DimensionedField<scalar, GeoMesh>> tAtan2 =
        reuseTmpTmpDimensionedField
        <
            scalar,
            scalar,
            scalar,
            GeoMesh,
            PrimitiveField1,
            PrimitiveField2
        >::New
        (
            tdsf1,
            tdsf2,
            "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
            atan2(dsf1.dimensions(), dsf2.dimensions())
        );

    atan2
    (
        tAtan2.ref().primitiveFieldRef(),
        dsf1.primitiveField(),
        dsf2.primitiveField()
    );

    tdsf1.clear();
    tdsf2.clear();

    return tAtan2;
}


template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf,
    const dimensionedScalar& ds
)
{
    tmp<DimensionedField<scalar, GeoMesh, Field>> tAtan2
    (
        DimensionedField<scalar, GeoMesh, Field>::New
        (
            "atan2(" + dsf.name() + ',' + ds.name() + ')',
            dsf.mesh(),
            atan2(dsf.dimensions(), ds)
        )
    );

    atan2(tAtan2.ref().primitiveFieldRef(), dsf.primitiveField(), ds.value());

    return tAtan2;
}

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField>>& tdsf,
    const dimensionedScalar& ds
)
{
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh, Field>> tAtan2 = New
    (
        tdsf,
        "atan2(" + dsf.name() + ',' + ds.name() + ')',
        atan2(dsf.dimensions(), ds)
    );

    atan2(tAtan2.ref().primitiveFieldRef(), dsf.primitiveField(), ds.value());

    tdsf.clear();

    return tAtan2;
}

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf,
    const scalar& s
)
{
    return atan2(dsf, dimensionedScalar(s));
}

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField>>& tdsf,
    const scalar& s
)
{
    return atan2(tdsf, dimensionedScalar(s));
}


template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const dimensionedScalar& ds,
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf
)
{
    tmp<DimensionedField<scalar, GeoMesh, Field>> tAtan2
    (
        DimensionedField<scalar, GeoMesh, Field>::New
        (
            "atan2(" + ds.name() + ',' + dsf.name() + ')',
            dsf.mesh(),
            atan2(ds, dsf.dimensions())
        )
    );

    atan2(tAtan2.ref().primitiveFieldRef(), ds.value(), dsf.primitiveField());

    return tAtan2;
}


template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const dimensionedScalar& ds,
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField>>& tdsf
)
{
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh, Field>> tAtan2 = New
    (
        tdsf,
        "atan2(" + ds.name() + ',' + dsf.name() + ')',
        atan2(ds, dsf.dimensions())
    );

    atan2(tAtan2.ref().primitiveFieldRef(), ds.value(), dsf.primitiveField());

    tdsf.clear();

    return tAtan2;
}

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const scalar& s,
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf
)
{
    return atan2(dimensionedScalar(s), dsf);
}

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> atan2
(
    const scalar& s,
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField>>& tdsf
)
{
    return atan2(dimensionedScalar(s), tdsf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, scalar, pow3, pow3)
UNARY_FUNCTION(scalar, scalar, pow4, pow4)
UNARY_FUNCTION(scalar, scalar, pow5, pow5)
UNARY_FUNCTION(scalar, scalar, pow6, pow6)
UNARY_FUNCTION(scalar, scalar, pow025, pow025)
UNARY_FUNCTION(scalar, scalar, sqrt, sqrt)
UNARY_FUNCTION(scalar, scalar, cbrt, cbrt)
UNARY_FUNCTION(scalar, scalar, sign, sign)
UNARY_FUNCTION(scalar, scalar, pos, pos)
UNARY_FUNCTION(scalar, scalar, pos0, pos0)
UNARY_FUNCTION(scalar, scalar, neg, neg)
UNARY_FUNCTION(scalar, scalar, neg0, neg0)
UNARY_FUNCTION(scalar, scalar, posPart, posPart)
UNARY_FUNCTION(scalar, scalar, negPart, negPart)

UNARY_FUNCTION(scalar, scalar, exp, trans)
UNARY_FUNCTION(scalar, scalar, log, trans)
UNARY_FUNCTION(scalar, scalar, log10, trans)
UNARY_FUNCTION(scalar, scalar, sin, trans)
UNARY_FUNCTION(scalar, scalar, cos, trans)
UNARY_FUNCTION(scalar, scalar, tan, trans)
UNARY_FUNCTION(scalar, scalar, asin, trans)
UNARY_FUNCTION(scalar, scalar, acos, trans)
UNARY_FUNCTION(scalar, scalar, atan, trans)
UNARY_FUNCTION(scalar, scalar, sinh, trans)
UNARY_FUNCTION(scalar, scalar, cosh, trans)
UNARY_FUNCTION(scalar, scalar, tanh, trans)
UNARY_FUNCTION(scalar, scalar, asinh, trans)
UNARY_FUNCTION(scalar, scalar, acosh, trans)
UNARY_FUNCTION(scalar, scalar, atanh, trans)
UNARY_FUNCTION(scalar, scalar, erf, trans)
UNARY_FUNCTION(scalar, scalar, erfc, trans)
UNARY_FUNCTION(scalar, scalar, lgamma, trans)
UNARY_FUNCTION(scalar, scalar, j0, trans)
UNARY_FUNCTION(scalar, scalar, j1, trans)
UNARY_FUNCTION(scalar, scalar, y0, trans)
UNARY_FUNCTION(scalar, scalar, y1, trans)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BesselFunc(func)                                                       \
                                                                               \
template<class GeoMesh, template<class> class PrimitiveField>                  \
tmp<DimensionedField<scalar, GeoMesh, Field>> func                             \
(                                                                              \
    const int n,                                                               \
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf               \
)                                                                              \
{                                                                              \
    if (!dsf.dimensions().dimensionless())                                     \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "dsf not dimensionless"                                         \
            << abort(FatalError);                                              \
    }                                                                          \
                                                                               \
    tmp<DimensionedField<scalar, GeoMesh, Field>> tFunc                        \
    (                                                                          \
        DimensionedField<scalar, GeoMesh, Field>::New                          \
        (                                                                      \
            #func "(" + name(n) + ',' + dsf.name() + ')',                      \
            dsf.mesh(),                                                        \
            dimless                                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    func(tFunc.ref().primitiveFieldRef(), n, dsf.primitiveField());            \
                                                                               \
    return tFunc;                                                              \
}                                                                              \
                                                                               \
template<class GeoMesh, template<class> class PrimitiveField>                  \
tmp<DimensionedField<scalar, GeoMesh, Field>> func                             \
(                                                                              \
    const int n,                                                               \
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField>>& tdsf         \
)                                                                              \
{                                                                              \
    const DimensionedField<scalar, GeoMesh, PrimitiveField>& dsf = tdsf();     \
                                                                               \
    if (!dsf.dimensions().dimensionless())                                     \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << " : dsf not dimensionless"                                      \
            << abort(FatalError);                                              \
    }                                                                          \
                                                                               \
    tmp<DimensionedField<scalar, GeoMesh, Field>> tFunc                        \
    (                                                                          \
        New                                                                    \
        (                                                                      \
            tdsf,                                                              \
            #func "(" + name(n) + ',' + dsf.name() + ')',                      \
            dimless                                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    func(tFunc.ref().primitiveFieldRef(), n, dsf.primitiveField());            \
                                                                               \
    tdsf.clear();                                                              \
                                                                               \
    return tFunc;                                                              \
}

BesselFunc(jn)
BesselFunc(yn)

#undef BesselFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
