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

#include "GeometricScalarField.H"

#define TEMPLATE                                                               \
    template                                                                   \
    <                                                                          \
        class GeoMesh,                                                         \
        template<class> class PrimitiveField                                   \
    >
#define TEMPLATE2                                                              \
    template                                                                   \
    <                                                                          \
        class GeoMesh,                                                         \
        template<class> class PrimitiveField1,                                 \
        template<class> class PrimitiveField2                                  \
    >
#define TEMPLATE3                                                              \
    template                                                                   \
    <                                                                          \
        class GeoMesh,                                                         \
        template<class> class PrimitiveField1,                                 \
        template<class> class PrimitiveField2,                                 \
        template<class> class PrimitiveField3                                  \
    >
#include "GeometricFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

TEMPLATE2
void stabilise
(
    GeometricField<scalar, GeoMesh, PrimitiveField1>& result,
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf,
    const dimensioned<scalar>& ds
)
{
    stabilise(result.primitiveFieldRef(), gsf.primitiveField(), ds.value());
    stabilise(result.boundaryFieldRef(), gsf.boundaryField(), ds.value());
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> stabilise
(
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf,
    const dimensioned<scalar>& ds
)
{
    tmp<GeometricField<scalar, GeoMesh, Field>> tRes
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "stabilise(" + gsf.name() + ',' + ds.name() + ')',
            gsf.mesh(),
            ds.dimensions() + gsf.dimensions()
        )
    );

    stabilise(tRes.ref(), gsf, ds);

    return tRes;
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> stabilise
(
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField>>& tgsf,
    const dimensioned<scalar>& ds
)
{
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf = tgsf();

    tmp<GeometricField<scalar, GeoMesh, Field>> tRes
    (
        New
        (
            tgsf,
            "stabilise(" + gsf.name() + ',' + ds.name() + ')',
            ds.dimensions() + gsf.dimensions()
        )
    );

    stabilise(tRes.ref(), gsf, ds);

    tgsf.clear();

    return tRes;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(scalar, scalar, scalar, +, '+', add)
BINARY_TYPE_OPERATOR(scalar, scalar, scalar, -, '-', subtract)

BINARY_OPERATOR(scalar, scalar, scalar, *, '*', multiply)
BINARY_OPERATOR(scalar, scalar, scalar, /, '|', divide)

BINARY_TYPE_OPERATOR_SF(scalar, scalar, scalar, /, '|', divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEMPLATE3
void pow
(
    GeometricField<scalar, GeoMesh, PrimitiveField1>& Pow,
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf1,
    const GeometricField<scalar, GeoMesh, PrimitiveField3>& gsf2
)
{
    pow(Pow.primitiveFieldRef(), gsf1.primitiveField(), gsf2.primitiveField());
    pow(Pow.boundaryFieldRef(), gsf1.boundaryField(), gsf2.boundaryField());
}


TEMPLATE2
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf1,
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf2
)
{
    if (!gsf1.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base field is not dimensionless: " << gsf1.dimensions()
            << exit(FatalError);
    }

    if (!gsf2.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << gsf2.dimensions()
            << exit(FatalError);
    }

    tmp<GeometricField<scalar, GeoMesh, Field>> tPow
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "pow(" + gsf1.name() + ',' + gsf2.name() + ')',
            gsf1.mesh(),
            dimless
        )
    );

    pow(tPow.ref(), gsf1, gsf2);

    return tPow;
}


TEMPLATE2
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField1>>& tgsf1,
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf2
)
{
    const GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf1 =
        tgsf1();

    if (!gsf1.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base field is not dimensionless: " << gsf1.dimensions()
            << exit(FatalError);
    }

    if (!gsf2.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << gsf2.dimensions()
            << exit(FatalError);
    }

    tmp<GeometricField<scalar, GeoMesh, Field>> tPow
    (
        reuseTmpGeometricField<scalar, scalar, GeoMesh, PrimitiveField1>::New
        (
            tgsf1,
            "pow(" + gsf1.name() + ',' + gsf2.name() + ')',
            dimless
        )
    );

    pow(tPow.ref(), gsf1, gsf2);

    tgsf1.clear();

    return tPow;
}


TEMPLATE2
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf1,
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField2>>& tgsf2
)
{
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf2 =
        tgsf2();

    if (!gsf1.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base field is not dimensionless: " << gsf1.dimensions()
            << exit(FatalError);
    }

    if (!gsf2.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << gsf2.dimensions()
            << exit(FatalError);
    }

    tmp<GeometricField<scalar, GeoMesh, Field>> tPow
    (
        reuseTmpGeometricField<scalar, scalar, GeoMesh, PrimitiveField2>::New
        (
            tgsf2,
            "pow(" + gsf1.name() + ',' + gsf2.name() + ')',
            dimless
        )
    );

    pow(tPow.ref(), gsf1, gsf2);

    tgsf2.clear();

    return tPow;
}


TEMPLATE2
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField1>>& tgsf1,
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField2>>& tgsf2
)
{
    const GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf1 = tgsf1();
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf2 = tgsf2();

    if (!gsf1.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base field is not dimensionless: " << gsf1.dimensions()
            << exit(FatalError);
    }

    if (!gsf2.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << gsf2.dimensions()
            << exit(FatalError);
    }

    tmp<GeometricField<scalar, GeoMesh, Field>> tPow
    (
        reuseTmpTmpGeometricField
        <
            scalar,
            scalar,
            scalar,
            GeoMesh,
            PrimitiveField1,
            PrimitiveField2
        >::New
        (
            tgsf1,
            tgsf2,
            "pow(" + gsf1.name() + ',' + gsf2.name() + ')',
            dimless
        )
    );

    pow(tPow.ref(), gsf1, gsf2);

    tgsf1.clear();
    tgsf2.clear();

    return tPow;
}


TEMPLATE2
void pow
(
    GeometricField<scalar, GeoMesh, PrimitiveField1>& tPow,
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf,
    const dimensioned<scalar>& ds
)
{
    pow(tPow.primitiveFieldRef(), gsf.primitiveField(), ds.value());
    pow(tPow.boundaryFieldRef(), gsf.boundaryField(), ds.value());
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf,
    const dimensionedScalar& ds
)
{
    if (!ds.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent is not dimensionless: " << ds.dimensions()
            << exit(FatalError);
    }

    tmp<GeometricField<scalar, GeoMesh, Field>> tPow
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "pow(" + gsf.name() + ',' + ds.name() + ')',
            gsf.mesh(),
            pow(gsf.dimensions(), ds)
        )
    );

    pow(tPow.ref(), gsf, ds);

    return tPow;
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField>>& tgsf,
    const dimensionedScalar& ds
)
{
    if (!ds.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent is not dimensionless: " << ds.dimensions()
            << exit(FatalError);
    }

    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf = tgsf();

    tmp<GeometricField<scalar, GeoMesh, Field>> tPow
    (
        reuseTmpGeometricField<scalar, scalar, GeoMesh, PrimitiveField>::New
        (
            tgsf,
            "pow(" + gsf.name() + ',' + ds.name() + ')',
            pow(gsf.dimensions(), ds)
        )
    );

    pow(tPow.ref(), gsf, ds);

    tgsf.clear();

    return tPow;
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf,
    const scalar& s
)
{
    return pow(gsf, dimensionedScalar(s));
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField>>& tgsf,
    const scalar& s
)
{
    return pow(tgsf, dimensionedScalar(s));
}


TEMPLATE2
void pow
(
    GeometricField<scalar, GeoMesh, PrimitiveField1>& tPow,
    const dimensioned<scalar>& ds,
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf
)
{
    pow(tPow.primitiveFieldRef(), ds.value(), gsf.primitiveField());
    pow(tPow.boundaryFieldRef(), ds.value(), gsf.boundaryField());
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const dimensionedScalar& ds,
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf
)
{
    if (!ds.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base scalar is not dimensionless: " << ds.dimensions()
            << exit(FatalError);
    }

    if (!gsf.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << gsf.dimensions()
            << exit(FatalError);
    }

    tmp<GeometricField<scalar, GeoMesh, Field>> tPow
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "pow(" + ds.name() + ',' + gsf.name() + ')',
            gsf.mesh(),
            dimless
        )
    );

    pow(tPow.ref(), ds, gsf);

    return tPow;
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const dimensionedScalar& ds,
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField>>& tgsf
)
{
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf = tgsf();

    if (!ds.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Base scalar is not dimensionless: " << ds.dimensions()
            << exit(FatalError);
    }

    if (!gsf.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent field is not dimensionless: " << gsf.dimensions()
            << exit(FatalError);
    }

    tmp<GeometricField<scalar, GeoMesh, Field>> tPow
    (
        reuseTmpGeometricField<scalar, scalar, GeoMesh, PrimitiveField>::New
        (
            tgsf,
            "pow(" + ds.name() + ',' + gsf.name() + ')',
            dimless
        )
    );

    pow(tPow.ref(), ds, gsf);

    tgsf.clear();

    return tPow;
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const scalar& s,
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf
)
{
    return pow(dimensionedScalar(s), gsf);
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> pow
(
    const scalar& s,
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField>>& tgsf
)
{
    return pow(dimensionedScalar(s), tgsf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEMPLATE3
void atan2
(
    GeometricField<scalar, GeoMesh, PrimitiveField1>& Atan2,
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf1,
    const GeometricField<scalar, GeoMesh, PrimitiveField3>& gsf2
)
{
    atan2
    (
        Atan2.primitiveFieldRef(),
        gsf1.primitiveField(),
        gsf2.primitiveField()
    );
    atan2
    (
        Atan2.boundaryFieldRef(),
        gsf1.boundaryField(),
        gsf2.boundaryField()
    );
}


TEMPLATE2
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf1,
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf2
)
{
    tmp<GeometricField<scalar, GeoMesh, Field>> tAtan2
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "atan2(" + gsf1.name() + ',' + gsf2.name() + ')',
            gsf1.mesh(),
            atan2(gsf1.dimensions(), gsf2.dimensions())
        )
    );

    atan2(tAtan2.ref(), gsf1, gsf2);

    return tAtan2;
}


TEMPLATE2
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField1>>& tgsf1,
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf2
)
{
    const GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf1 = tgsf1();

    tmp<GeometricField<scalar, GeoMesh, Field>> tAtan2
    (
        reuseTmpGeometricField
        <scalar, scalar, GeoMesh, PrimitiveField1>::New
        (
            tgsf1,
            "atan2(" + gsf1.name() + ',' + gsf2.name() + ')',
            atan2(gsf1.dimensions(), gsf2.dimensions())
        )
    );

    atan2(tAtan2.ref(), gsf1, gsf2);

    tgsf1.clear();

    return tAtan2;
}


TEMPLATE2
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf1,
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField2>>& tgsf2
)
{
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf2 = tgsf2();

    tmp<GeometricField<scalar, GeoMesh, Field>> tAtan2
    (
        reuseTmpGeometricField
        <scalar, scalar, GeoMesh, PrimitiveField2>::New
        (
            tgsf2,
            "atan2(" + gsf1.name() + ',' + gsf2.name() + ')',
            atan2( gsf1.dimensions(), gsf2.dimensions())
        )
    );

    atan2(tAtan2.ref(), gsf1, gsf2);

    tgsf2.clear();

    return tAtan2;
}


TEMPLATE2
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField1>>& tgsf1,
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField2>>& tgsf2
)
{
    const GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf1 = tgsf1();
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf2 = tgsf2();

    tmp<GeometricField<scalar, GeoMesh, Field>> tAtan2
    (
        reuseTmpTmpGeometricField
        <
            scalar,
            scalar,
            scalar,
            GeoMesh,
            PrimitiveField1,
            PrimitiveField2
        >::New
        (
            tgsf1,
            tgsf2,
            "atan2(" + gsf1.name() + ',' + gsf2.name() + ')',
            atan2(gsf1.dimensions(), gsf2.dimensions())
        )
    );

    atan2(tAtan2.ref(), gsf1, gsf2);

    tgsf1.clear();
    tgsf2.clear();

    return tAtan2;
}


TEMPLATE2
void atan2
(
    GeometricField<scalar, GeoMesh, PrimitiveField1>& tAtan2,
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf,
    const dimensioned<scalar>& ds
)
{
    atan2(tAtan2.primitiveFieldRef(), gsf.primitiveField(), ds.value());
    atan2(tAtan2.boundaryFieldRef(), gsf.boundaryField(), ds.value());
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf,
    const dimensionedScalar& ds
)
{
    tmp<GeometricField<scalar, GeoMesh, Field>> tAtan2
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "atan2(" + gsf.name() + ',' + ds.name() + ')',
            gsf.mesh(),
            atan2(gsf.dimensions(), ds)
        )
    );

    atan2(tAtan2.ref(), gsf, ds);

    return tAtan2;
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField>>& tgsf,
    const dimensionedScalar& ds
)
{
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf = tgsf();

    tmp<GeometricField<scalar, GeoMesh, Field>> tAtan2
    (
        reuseTmpGeometricField<scalar, scalar, GeoMesh, PrimitiveField>::New
        (
            tgsf,
            "atan2(" + gsf.name() + ',' + ds.name() + ')',
            atan2(gsf.dimensions(), ds)
        )
    );

    atan2(tAtan2.ref(), gsf, ds);

    tgsf.clear();

    return tAtan2;
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf,
    const scalar& s
)
{
    return atan2(gsf, dimensionedScalar(s));
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField>>& tgsf,
    const scalar& s
)
{
    return atan2(tgsf, dimensionedScalar(s));
}


TEMPLATE2
void atan2
(
    GeometricField<scalar, GeoMesh, PrimitiveField1>& tAtan2,
    const dimensioned<scalar>& ds,
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf
)
{
    atan2(tAtan2.primitiveFieldRef(), ds.value(), gsf.primitiveField());
    atan2(tAtan2.boundaryFieldRef(), ds.value(), gsf.boundaryField());
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const dimensionedScalar& ds,
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf
)
{
    tmp<GeometricField<scalar, GeoMesh, Field>> tAtan2
    (
        GeometricField<scalar, GeoMesh, Field>::New
        (
            "atan2(" + ds.name() + ',' + gsf.name() + ')',
            gsf.mesh(),
            atan2(ds, gsf.dimensions())
        )
    );

    atan2(tAtan2.ref(), ds, gsf);

    return tAtan2;
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const dimensionedScalar& ds,
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField>>& tgsf
)
{
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf =tgsf();

    tmp<GeometricField<scalar, GeoMesh, Field>> tAtan2
    (
        reuseTmpGeometricField<scalar, scalar, GeoMesh, PrimitiveField>::New
        (
            tgsf,
            "atan2(" + ds.name() + ',' + gsf.name() + ')',
            atan2(ds, gsf.dimensions())
        )
    );

    atan2(tAtan2.ref(), ds, gsf);

    tgsf.clear();

    return tAtan2;
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const scalar& s,
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf
)
{
    return atan2(dimensionedScalar(s), gsf);
}


TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> atan2
(
    const scalar& s,
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField>>& tgsf
)
{
    return atan2(dimensionedScalar(s), tgsf);
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
TEMPLATE2                                                                      \
void func                                                                      \
(                                                                              \
    GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf,                     \
    const int n,                                                               \
    const GeometricField<scalar, GeoMesh, PrimitiveField2>& gsf1               \
)                                                                              \
{                                                                              \
    func(gsf.primitiveFieldRef(), n, gsf1.primitiveField());                   \
    func(gsf.boundaryFieldRef(), n, gsf1.boundaryField());                     \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<scalar, GeoMesh, Field>> func                               \
(                                                                              \
    const int n,                                                               \
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf                 \
)                                                                              \
{                                                                              \
    if (!gsf.dimensions().dimensionless())                                     \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "gsf not dimensionless"                                         \
            << abort(FatalError);                                              \
    }                                                                          \
                                                                               \
    tmp<GeometricField<scalar, GeoMesh, Field>> tFunc                          \
    (                                                                          \
        GeometricField<scalar, GeoMesh, Field>::New                            \
        (                                                                      \
            #func "(" + gsf.name() + ')',                                      \
            gsf.mesh(),                                                        \
            dimless                                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    func(tFunc.ref(), n, gsf);                                                 \
                                                                               \
    return tFunc;                                                              \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<scalar, GeoMesh, Field>> func                               \
(                                                                              \
    const int n,                                                               \
    const tmp<GeometricField<scalar, GeoMesh, PrimitiveField>>& tgsf           \
)                                                                              \
{                                                                              \
    const GeometricField<scalar, GeoMesh, PrimitiveField>& gsf = tgsf();       \
                                                                               \
    if (!gsf.dimensions().dimensionless())                                     \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << " : gsf not dimensionless"                                      \
            << abort(FatalError);                                              \
    }                                                                          \
                                                                               \
    tmp<GeometricField<scalar, GeoMesh, Field>> tFunc                          \
    (                                                                          \
        reuseTmpGeometricField<scalar, scalar, GeoMesh, PrimitiveField>::New   \
        (                                                                      \
            tgsf,                                                              \
            #func "(" + gsf.name() + ')',                                      \
            dimless                                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    func(tFunc.ref(), n, gsf);                                                 \
                                                                               \
    tgsf.clear();                                                              \
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
