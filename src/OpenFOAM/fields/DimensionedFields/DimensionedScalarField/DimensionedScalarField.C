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

#include "DimensionedScalarField.H"

#define TEMPLATE template<class GeoMesh>
#include "DimensionedFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > stabilise
(
    const DimensionedField<scalar, GeoMesh>& dsf,
    const dimensioned<scalar>& ds
)
{
    tmp<DimensionedField<scalar, GeoMesh> > tRes
    (
        new DimensionedField<scalar, GeoMesh>
        (
            IOobject
            (
                "stabilise(" + dsf.name() + ',' + ds.name() + ')',
                dsf.instance(),
                dsf.db()
            ),
            dsf.mesh(),
            dsf.dimensions() + ds.dimensions()
        )
    );

    stabilise(tRes().field(), dsf.field(), ds.value());

    return tRes;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > stabilise
(
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf,
    const dimensioned<scalar>& ds
)
{
    const DimensionedField<scalar, GeoMesh>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh> > tRes =
        reuseTmpDimensionedField<scalar, scalar, GeoMesh>::New
        (
            tdsf,
            "stabilise(" + dsf.name() + ',' + ds.name() + ')',
            dsf.dimensions() + ds.dimensions()
        );

    stabilise(tRes().field(), dsf.field(), ds.value());

    reuseTmpDimensionedField<scalar, scalar, GeoMesh>::clear(tdsf);

    return tRes;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(scalar, scalar, scalar, +, '+', add)
BINARY_TYPE_OPERATOR(scalar, scalar, scalar, -, '-', subtract)

BINARY_OPERATOR(scalar, scalar, scalar, *, '*', multiply)
BINARY_OPERATOR(scalar, scalar, scalar, /, '|', divide)

BINARY_TYPE_OPERATOR_SF(scalar, scalar, scalar, /, '|', divide)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const DimensionedField<scalar, GeoMesh>& dsf1,
    const DimensionedField<scalar, GeoMesh>& dsf2
)
{
    tmp<DimensionedField<scalar, GeoMesh> > tPow
    (
        new DimensionedField<scalar, GeoMesh>
        (
            IOobject
            (
                "pow(" + dsf1.name() + ',' + dsf2.name() + ')',
                dsf1.instance(),
                dsf1.db()
            ),
            dsf1.mesh(),
            pow
            (
                dsf1.dimensions(),
                dimensionedScalar("1", 1.0, dsf2.dimensions())
            )
        )
    );

    pow(tPow().field(), dsf1.field(), dsf2.field());

    return tPow;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf1,
    const DimensionedField<scalar, GeoMesh>& dsf2
)
{
    const DimensionedField<scalar, GeoMesh>& dsf1 = tdsf1();

    tmp<DimensionedField<scalar, GeoMesh> > tPow =
        reuseTmpDimensionedField<scalar, scalar, GeoMesh>::New
        (
            tdsf1,
            "pow(" + dsf1.name() + ',' + dsf2.name() + ')',
            pow
            (
                dsf1.dimensions(),
                dimensionedScalar("1", 1.0, dsf2.dimensions())
            )
        );

    pow(tPow().field(), dsf1.field(), dsf2.field());

    reuseTmpDimensionedField<scalar, scalar, GeoMesh>::clear(tdsf1);

    return tPow;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const DimensionedField<scalar, GeoMesh>& dsf1,
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf2
)
{
    const DimensionedField<scalar, GeoMesh>& dsf2 = tdsf2();

    tmp<DimensionedField<scalar, GeoMesh> > tPow =
        reuseTmpDimensionedField<scalar, scalar, GeoMesh>::New
        (
            tdsf2,
            "pow(" + dsf1.name() + ',' + dsf2.name() + ')',
            pow
            (
                dsf1.dimensions(),
                dimensionedScalar("1", 1.0, dsf2.dimensions())
            )
        );

    pow(tPow().field(), dsf1.field(), dsf2.field());

    reuseTmpDimensionedField<scalar, scalar, GeoMesh>::clear(tdsf2);

    return tPow;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf1,
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf2
)
{
    const DimensionedField<scalar, GeoMesh>& dsf1 = tdsf1();
    const DimensionedField<scalar, GeoMesh>& dsf2 = tdsf2();

    tmp<DimensionedField<scalar, GeoMesh> > tPow =
        reuseTmpTmpDimensionedField<scalar, scalar, scalar, scalar, GeoMesh>::
        New
        (
            tdsf1,
            tdsf2,
            "pow(" + dsf1.name() + ',' + dsf2.name() + ')',
            pow
            (
                dsf1.dimensions(),
                dimensionedScalar("1", 1.0, dsf2.dimensions())
            )
        );

    pow(tPow().field(), dsf1.field(), dsf2.field());

    reuseTmpTmpDimensionedField<scalar, scalar, scalar, scalar, GeoMesh>::clear
    (
        tdsf1,
        tdsf2
    );

    return tPow;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const DimensionedField<scalar, GeoMesh>& dsf,
    const dimensionedScalar& ds
)
{
    tmp<DimensionedField<scalar, GeoMesh> > tPow
    (
        new DimensionedField<scalar, GeoMesh>
        (
            IOobject
            (
                "pow(" + dsf.name() + ',' + ds.name() + ')',
                dsf.instance(),
                dsf.db()
            ),
            dsf.mesh(),
            pow(dsf.dimensions(), ds)
        )
    );

    pow(tPow().field(), dsf.field(), ds.value());

    return tPow;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf,
    const dimensionedScalar& ds
)
{
    const DimensionedField<scalar, GeoMesh>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh> > tPow =
        reuseTmpDimensionedField<scalar, scalar, GeoMesh>::New
        (
            tdsf,
            "pow(" + dsf.name() + ',' + ds.name() + ')',
            pow(dsf.dimensions(), ds)
        );

    pow(tPow().field(), dsf.field(), ds.value());

    reuseTmpDimensionedField<scalar, scalar, GeoMesh>::clear(tdsf);

    return tPow;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const DimensionedField<scalar, GeoMesh>& dsf,
    const scalar& s
)
{
    return pow(dsf, dimensionedScalar(s));
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf,
    const scalar& s
)
{
    return pow(tdsf, dimensionedScalar(s));
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const dimensionedScalar& ds,
    const DimensionedField<scalar, GeoMesh>& dsf
)
{
    tmp<DimensionedField<scalar, GeoMesh> > tPow
    (
        new DimensionedField<scalar, GeoMesh>
        (
            IOobject
            (
                "pow(" + ds.name() + ',' + dsf.name() + ')',
                dsf.instance(),
                dsf.db()
            ),
            dsf.mesh(),
            pow(ds, dsf.dimensions())
        )
    );

    pow(tPow().field(), ds.value(), dsf.field());

    return tPow;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const dimensionedScalar& ds,
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf
)
{
    const DimensionedField<scalar, GeoMesh>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh> > tPow =
        reuseTmpDimensionedField<scalar, scalar, GeoMesh>::New
        (
            tdsf,
            "pow(" + ds.name() + ',' + dsf.name() + ')',
            pow(ds, dsf.dimensions())
        );

    pow(tPow().field(), ds.value(), dsf.field());

    reuseTmpDimensionedField<scalar, scalar, GeoMesh>::clear(tdsf);

    return tPow;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const scalar& s,
    const DimensionedField<scalar, GeoMesh>& dsf
)
{
    return pow(dimensionedScalar(s), dsf);
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > pow
(
    const scalar& s,
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf
)
{
    return pow(dimensionedScalar(s), tdsf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const DimensionedField<scalar, GeoMesh>& dsf1,
    const DimensionedField<scalar, GeoMesh>& dsf2
)
{
    tmp<DimensionedField<scalar, GeoMesh> > tAtan2
    (
        new DimensionedField<scalar, GeoMesh>
        (
            IOobject
            (
                "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
                dsf1.instance(),
                dsf1.db()
            ),
            dsf1.mesh(),
            atan2(dsf1.dimensions(), dsf2.dimensions())
        )
    );

    atan2(tAtan2().field(), dsf1.field(), dsf2.field());

    return tAtan2;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf1,
    const DimensionedField<scalar, GeoMesh>& dsf2
)
{
    const DimensionedField<scalar, GeoMesh>& dsf1 = tdsf1();

    tmp<DimensionedField<scalar, GeoMesh> > tAtan2 =
        reuseTmpDimensionedField<scalar, scalar, GeoMesh>::New
        (
            tdsf1,
            "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
            atan2(dsf1.dimensions(), dsf2.dimensions())
        );

    atan2(tAtan2().field(), dsf1.field(), dsf2.field());

    reuseTmpDimensionedField<scalar, scalar, GeoMesh>::clear(tdsf1);

    return tAtan2;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const DimensionedField<scalar, GeoMesh>& dsf1,
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf2
)
{
    const DimensionedField<scalar, GeoMesh>& dsf2 = tdsf2();

    tmp<DimensionedField<scalar, GeoMesh> > tAtan2 =
        reuseTmpDimensionedField<scalar, scalar, GeoMesh>::New
        (
            tdsf2,
            "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
            atan2(dsf1.dimensions(), dsf2.dimensions())
        );

    atan2(tAtan2().field(), dsf1.field(), dsf2.field());

    reuseTmpDimensionedField<scalar, scalar, GeoMesh>::clear(tdsf2);

    return tAtan2;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf1,
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf2
)
{
    const DimensionedField<scalar, GeoMesh>& dsf1 = tdsf1();
    const DimensionedField<scalar, GeoMesh>& dsf2 = tdsf2();

    tmp<DimensionedField<scalar, GeoMesh> > tAtan2 =
        reuseTmpTmpDimensionedField<scalar, scalar, scalar, scalar, GeoMesh>::
        New
        (
            tdsf1,
            tdsf2,
            "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
            atan2(dsf1.dimensions(), dsf2.dimensions())
        );

    atan2(tAtan2().field(), dsf1.field(), dsf2.field());

    reuseTmpTmpDimensionedField<scalar, scalar, scalar, scalar, GeoMesh>::clear
    (
        tdsf1,
        tdsf2
    );

    return tAtan2;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const DimensionedField<scalar, GeoMesh>& dsf,
    const dimensionedScalar& ds
)
{
    tmp<DimensionedField<scalar, GeoMesh> > tAtan2
    (
        new DimensionedField<scalar, GeoMesh>
        (
            IOobject
            (
                "atan2(" + dsf.name() + ',' + ds.name() + ')',
                dsf.instance(),
                dsf.db()
            ),
            dsf.mesh(),
            atan2(dsf.dimensions(), ds)
        )
    );

    atan2(tAtan2().field(), dsf.field(), ds.value());

    return tAtan2;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf,
    const dimensionedScalar& ds
)
{
    const DimensionedField<scalar, GeoMesh>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh> > tAtan2 =
        reuseTmpDimensionedField<scalar, scalar, GeoMesh>::New
        (
            tdsf,
            "atan2(" + dsf.name() + ',' + ds.name() + ')',
            atan2(dsf.dimensions(), ds)
        );

    atan2(tAtan2().field(), dsf.field(), ds.value());

    reuseTmpDimensionedField<scalar, scalar, GeoMesh>::clear(tdsf);

    return tAtan2;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const DimensionedField<scalar, GeoMesh>& dsf,
    const scalar& s
)
{
    return atan2(dsf, dimensionedScalar(s));
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf,
    const scalar& s
)
{
    return atan2(tdsf, dimensionedScalar(s));
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const dimensionedScalar& ds,
    const DimensionedField<scalar, GeoMesh>& dsf
)
{
    tmp<DimensionedField<scalar, GeoMesh> > tAtan2
    (
        new DimensionedField<scalar, GeoMesh>
        (
            IOobject
            (
                "atan2(" + ds.name() + ',' + dsf.name() + ')',
                dsf.instance(),
                dsf.db()
            ),
            dsf.mesh(),
            atan2(ds, dsf.dimensions())
        )
    );

    atan2(tAtan2().field(), ds.value(), dsf.field());

    return tAtan2;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const dimensionedScalar& ds,
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf
)
{
    const DimensionedField<scalar, GeoMesh>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh> > tAtan2 =
        reuseTmpDimensionedField<scalar, scalar, GeoMesh>::New
        (
            tdsf,
            "atan2(" + ds.name() + ',' + dsf.name() + ')',
            atan2(ds, dsf.dimensions())
        );

    atan2(tAtan2().field(), ds.value(), dsf.field());

    reuseTmpDimensionedField<scalar, scalar, GeoMesh>::clear(tdsf);

    return tAtan2;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const scalar& s,
    const DimensionedField<scalar, GeoMesh>& dsf
)
{
    return atan2(dimensionedScalar(s), dsf);
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > atan2
(
    const scalar& s,
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf
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

#define BesselFunc(func)                                                    \
                                                                            \
template<class GeoMesh>                                                     \
tmp<DimensionedField<scalar, GeoMesh> > func                                \
(                                                                           \
    const int n,                                                            \
    const DimensionedField<scalar, GeoMesh>& dsf                            \
)                                                                           \
{                                                                           \
    if (!dsf.dimensions().dimensionless())                                  \
    {                                                                       \
        FatalErrorIn                                                        \
        (                                                                   \
            #func"(const int n, "                                           \
            "const DimensionedField<scalar, GeoMesh>& dsf)"                 \
        )   << "dsf not dimensionless"                                      \
            << abort(FatalError);                                           \
    }                                                                       \
                                                                            \
    tmp<DimensionedField<scalar, GeoMesh> > tFunc                           \
    (                                                                       \
        new DimensionedField<scalar, GeoMesh>                               \
        (                                                                   \
            IOobject                                                        \
            (                                                               \
                #func "(" + name(n) + ',' + dsf.name() + ')',               \
                dsf.instance(),                                             \
                dsf.db()                                                    \
            ),                                                              \
            dsf.mesh(),                                                     \
            dimless                                                         \
        )                                                                   \
    );                                                                      \
                                                                            \
    func(tFunc().field(), n, dsf.field());                                  \
                                                                            \
    return tFunc;                                                           \
}                                                                           \
                                                                            \
template<class GeoMesh>                                                     \
tmp<DimensionedField<scalar, GeoMesh> > func                                \
(                                                                           \
    const int n,                                                            \
    const tmp<DimensionedField<scalar, GeoMesh> >& tdsf                     \
)                                                                           \
{                                                                           \
    const DimensionedField<scalar, GeoMesh>& dsf = tdsf();                  \
                                                                            \
    if (!dsf.dimensions().dimensionless())                                  \
    {                                                                       \
        FatalErrorIn                                                        \
        (                                                                   \
            #func"(const int n, "                                           \
            "const tmp<DimensionedField<scalar, GeoMesh> >& dsf)"           \
        )   << " : dsf not dimensionless"                                   \
            << abort(FatalError);                                           \
    }                                                                       \
                                                                            \
    tmp<DimensionedField<scalar, GeoMesh> > tFunc                           \
    (                                                                       \
        reuseTmpDimensionedField<scalar, scalar, GeoMesh>::New              \
        (                                                                   \
            tdsf,                                                           \
            #func "(" + name(n) + ',' + dsf.name() + ')',                   \
            dimless                                                         \
        )                                                                   \
    );                                                                      \
                                                                            \
    func(tFunc().field(), n, dsf.field());                                  \
                                                                            \
    reuseTmpDimensionedField<scalar, scalar, GeoMesh>::clear(tdsf);         \
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
