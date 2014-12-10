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

#include "GeometricScalarField.H"

#define TEMPLATE template<template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<template<class> class PatchField, class GeoMesh>
void stabilise
(
    GeometricField<scalar, PatchField, GeoMesh>& result,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const dimensioned<scalar>& ds
)
{
    stabilise(result.internalField(), gsf.internalField(), ds.value());
    stabilise(result.boundaryField(), gsf.boundaryField(), ds.value());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > stabilise
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const dimensioned<scalar>& ds
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tRes
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "stabilise(" + gsf.name() + ',' + ds.name() + ')',
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf.mesh(),
            ds.dimensions() + gsf.dimensions()
        )
    );

    stabilise(tRes(), gsf, ds);

    return tRes;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > stabilise
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf,
    const dimensioned<scalar>& ds
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tRes
    (
        reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::New
        (
            tgsf,
            "stabilise(" + gsf.name() + ',' + ds.name() + ')',
            ds.dimensions() + gsf.dimensions()
        )
    );

    stabilise(tRes(), gsf, ds);

    reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::clear(tgsf);

    return tRes;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(scalar, scalar, scalar, +, '+', add)
BINARY_TYPE_OPERATOR(scalar, scalar, scalar, -, '-', subtract)

BINARY_OPERATOR(scalar, scalar, scalar, *, '*', multiply)
BINARY_OPERATOR(scalar, scalar, scalar, /, '|', divide)

BINARY_TYPE_OPERATOR_SF(scalar, scalar, scalar, /, '|', divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class PatchField, class GeoMesh>
void pow
(
    GeometricField<scalar, PatchField, GeoMesh>& Pow,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2
)
{
    pow(Pow.internalField(), gsf1.internalField(), gsf2.internalField());
    pow(Pow.boundaryField(), gsf1.boundaryField(), gsf2.boundaryField());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "pow(" + gsf1.name() + ',' + gsf2.name() + ')',
                gsf1.instance(),
                gsf1.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf1.mesh(),
            pow
            (
                gsf1.dimensions(),
                dimensionedScalar("1", gsf2.dimensions(), 1.0)
            )
        )
    );

    pow(tPow(), gsf1, gsf2);

    return tPow;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf1,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1 = tgsf1();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::New
        (
            tgsf1,
            "pow(" + gsf1.name() + ',' + gsf2.name() + ')',
            pow
            (
                gsf1.dimensions(),
                dimensionedScalar("1", gsf2.dimensions(), 1.0)
            )
        )
    );

    pow(tPow(), gsf1, gsf2);

    reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::clear(tgsf1);

    return tPow;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf2
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2 = tgsf2();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::New
        (
            tgsf2,
            "pow(" + gsf1.name() + ',' + gsf2.name() + ')',
            pow
            (
                gsf1.dimensions(),
                dimensionedScalar("1", gsf2.dimensions(), 1.0)
            )
        )
    );

    pow(tPow(), gsf1, gsf2);

    reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::clear(tgsf2);

    return tPow;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf1,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf2
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1 = tgsf1();
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2 = tgsf2();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        reuseTmpTmpGeometricField
            <scalar, scalar, scalar, scalar, PatchField, GeoMesh>::New
        (
            tgsf1,
            tgsf2,
            "pow(" + gsf1.name() + ',' + gsf2.name() + ')',
            pow
            (
                gsf1.dimensions(),
                dimensionedScalar("1", gsf2.dimensions(), 1.0)
            )
        )
    );

    pow(tPow(), gsf1, gsf2);

    reuseTmpTmpGeometricField
        <scalar, scalar, scalar, scalar, PatchField, GeoMesh>
        ::clear(tgsf1, tgsf2);

    return tPow;
}


template<template<class> class PatchField, class GeoMesh>
void pow
(
    GeometricField<scalar, PatchField, GeoMesh>& tPow,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const dimensioned<scalar>& ds
)
{
    pow(tPow.internalField(), gsf.internalField(), ds.value());
    pow(tPow.boundaryField(), gsf.boundaryField(), ds.value());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const dimensionedScalar& ds
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "pow(" + gsf.name() + ',' + ds.name() + ')',
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf.mesh(),
            pow(gsf.dimensions(), ds)
        )
    );

    pow(tPow(), gsf, ds);

    return tPow;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf,
    const dimensionedScalar& ds
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::New
        (
            tgsf,
            "pow(" + gsf.name() + ',' + ds.name() + ')',
            pow(gsf.dimensions(), ds)
        )
    );

    pow(tPow(), gsf, ds);

    reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::clear(tgsf);

    return tPow;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const scalar& s
)
{
    return pow(gsf, dimensionedScalar(s));
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf,
    const scalar& s
)
{
    return pow(tgsf, dimensionedScalar(s));
}


template<template<class> class PatchField, class GeoMesh>
void pow
(
    GeometricField<scalar, PatchField, GeoMesh>& tPow,
    const dimensioned<scalar>& ds,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    pow(tPow.internalField(), ds.value(), gsf.internalField());
    pow(tPow.boundaryField(), ds.value(), gsf.boundaryField());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const dimensionedScalar& ds,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "pow(" + ds.name() + ',' + gsf.name() + ')',
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf.mesh(),
            pow(ds, gsf.dimensions())
        )
    );

    pow(tPow(), ds, gsf);

    return tPow;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const dimensionedScalar& ds,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::New
        (
            tgsf,
            "pow(" + ds.name() + ',' + gsf.name() + ')',
            pow(ds, gsf.dimensions())
        )
    );

    pow(tPow(), ds, gsf);

    reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::clear(tgsf);

    return tPow;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const scalar& s,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    return pow(dimensionedScalar(s), gsf);
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const scalar& s,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf
)
{
    return pow(dimensionedScalar(s), tgsf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class PatchField, class GeoMesh>
void atan2
(
    GeometricField<scalar, PatchField, GeoMesh>& Atan2,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2
)
{
    atan2(Atan2.internalField(), gsf1.internalField(), gsf2.internalField());
    atan2(Atan2.boundaryField(), gsf1.boundaryField(), gsf2.boundaryField());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tAtan2
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "atan2(" + gsf1.name() + ',' + gsf2.name() + ')',
                gsf1.instance(),
                gsf1.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf1.mesh(),
            atan2(gsf1.dimensions(), gsf2.dimensions())
        )
    );

    atan2(tAtan2(), gsf1, gsf2);

    return tAtan2;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf1,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1 = tgsf1();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tAtan2
    (
        reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::New
        (
            tgsf1,
            "atan2(" + gsf1.name() + ',' + gsf2.name() + ')',
            atan2(gsf1.dimensions(), gsf2.dimensions())
        )
    );

    atan2(tAtan2(), gsf1, gsf2);

    reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::clear(tgsf1);

    return tAtan2;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf2
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2 = tgsf2();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tAtan2
    (
        reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::New
        (
            tgsf2,
            "atan2(" + gsf1.name() + ',' + gsf2.name() + ')',
            atan2( gsf1.dimensions(), gsf2.dimensions())
        )
    );

    atan2(tAtan2(), gsf1, gsf2);

    reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::clear(tgsf2);

    return tAtan2;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf1,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf2
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1 = tgsf1();
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2 = tgsf2();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tAtan2
    (
        reuseTmpTmpGeometricField
            <scalar, scalar, scalar, scalar, PatchField, GeoMesh>::New
        (
            tgsf1,
            tgsf2,
            "atan2(" + gsf1.name() + ',' + gsf2.name() + ')',
            atan2(gsf1.dimensions(), gsf2.dimensions())
        )
    );

    atan2(tAtan2(), gsf1, gsf2);

    reuseTmpTmpGeometricField
        <scalar, scalar, scalar, scalar, PatchField, GeoMesh>
        ::clear(tgsf1, tgsf2);

    return tAtan2;
}


template<template<class> class PatchField, class GeoMesh>
void atan2
(
    GeometricField<scalar, PatchField, GeoMesh>& tAtan2,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const dimensioned<scalar>& ds
)
{
    atan2(tAtan2.internalField(), gsf.internalField(), ds.value());
    atan2(tAtan2.boundaryField(), gsf.boundaryField(), ds.value());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const dimensionedScalar& ds
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tAtan2
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "atan2(" + gsf.name() + ',' + ds.name() + ')',
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf.mesh(),
            atan2(gsf.dimensions(), ds)
        )
    );

    atan2(tAtan2(), gsf, ds);

    return tAtan2;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf,
    const dimensionedScalar& ds
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tAtan2
    (
        reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::New
        (
            tgsf,
            "atan2(" + gsf.name() + ',' + ds.name() + ')',
            atan2(gsf.dimensions(), ds)
        )
    );

    atan2(tAtan2(), gsf, ds);

    reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::clear(tgsf);

    return tAtan2;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const scalar& s
)
{
    return atan2(gsf, dimensionedScalar(s));
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf,
    const scalar& s
)
{
    return atan2(tgsf, dimensionedScalar(s));
}


template<template<class> class PatchField, class GeoMesh>
void atan2
(
    GeometricField<scalar, PatchField, GeoMesh>& tAtan2,
    const dimensioned<scalar>& ds,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    atan2(tAtan2.internalField(), ds.value(), gsf.internalField());
    atan2(tAtan2.boundaryField(), ds.value(), gsf.boundaryField());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const dimensionedScalar& ds,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tAtan2
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "atan2(" + ds.name() + ',' + gsf.name() + ')',
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf.mesh(),
            atan2(ds, gsf.dimensions())
        )
    );

    atan2(tAtan2(), ds, gsf);

    return tAtan2;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const dimensionedScalar& ds,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tAtan2
    (
        reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::New
        (
            tgsf,
            "atan2(" + ds.name() + ',' + gsf.name() + ')',
            atan2(ds, gsf.dimensions())
        )
    );

    atan2(tAtan2(), ds, gsf);

    reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::clear(tgsf);

    return tAtan2;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const scalar& s,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    return atan2(dimensionedScalar(s), gsf);
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > atan2
(
    const scalar& s,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf
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
UNARY_FUNCTION(scalar, scalar, neg, neg)

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

#define BesselFunc(func)                                                     \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
void func                                                                   \
(                                                                           \
    GeometricField<scalar, PatchField, GeoMesh>& gsf,                       \
    const int n,                                                            \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1                 \
)                                                                           \
{                                                                           \
    func(gsf.internalField(), n, gsf1.internalField());                     \
    func(gsf.boundaryField(), n, gsf1.boundaryField());                     \
}                                                                           \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                      \
(                                                                           \
    const int n,                                                            \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf                  \
)                                                                           \
{                                                                           \
    if (!gsf.dimensions().dimensionless())                                  \
    {                                                                       \
        FatalErrorIn                                                        \
        (                                                                   \
            #func"(const int n, "                                           \
            "const GeometricField<scalar, PatchField, GeoMesh>& gsf)"       \
        )   << "gsf not dimensionless"                                      \
            << abort(FatalError);                                           \
    }                                                                       \
                                                                            \
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tFunc                 \
    (                                                                       \
        new GeometricField<scalar, PatchField, GeoMesh>                     \
        (                                                                   \
            IOobject                                                        \
            (                                                               \
                #func "(" + gsf.name() + ')',                               \
                gsf.instance(),                                             \
                gsf.db(),                                                   \
                IOobject::NO_READ,                                          \
                IOobject::NO_WRITE                                          \
            ),                                                              \
            gsf.mesh(),                                                     \
            dimless                                                         \
        )                                                                   \
    );                                                                      \
                                                                            \
    func(tFunc(), n, gsf);                                                  \
                                                                            \
    return tFunc;                                                           \
}                                                                           \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                      \
(                                                                           \
    const int n,                                                            \
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf           \
)                                                                           \
{                                                                           \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();        \
                                                                            \
    if (!gsf.dimensions().dimensionless())                                  \
    {                                                                       \
        FatalErrorIn                                                        \
        (                                                                   \
            #func"(const int n, "                                           \
            "const tmp<GeometricField<scalar, PatchField, GeoMesh> >& gsf)" \
        )   << " : gsf not dimensionless"                                   \
            << abort(FatalError);                                           \
    }                                                                       \
                                                                            \
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tFunc                 \
    (                                                                       \
        reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>::New    \
        (                                                                   \
            tgsf,                                                           \
            #func "(" + gsf.name() + ')',                                   \
            dimless                                                         \
        )                                                                   \
    );                                                                      \
                                                                            \
    func(tFunc(), n, gsf);                                                  \
                                                                            \
    reuseTmpGeometricField<scalar, scalar, PatchField, GeoMesh>             \
    ::clear(tgsf);                                                          \
                                                                            \
    return tFunc;                                                           \
}

BesselFunc(jn)
BesselFunc(yn)

#undef BesselFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
