/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#define TEMPLATE template<class Type, class GeoMesh>
#include "DimensionedFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, int r>
tmp<DimensionedField<typename powProduct<Type, r>::type, GeoMesh> >
pow
(
    const DimensionedField<Type, GeoMesh>& df,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;

    tmp<DimensionedField<powProductType, GeoMesh> > tPow
    (
        new DimensionedField<powProductType, GeoMesh>
        (
            IOobject
            (
                "pow(" + df.name() + ',' + name(r) + ')',
                df.instance(),
                df.db()
            ),
            df.mesh(),
            pow(df.dimensions(), r)
        )
    );

    pow<Type, r, GeoMesh>(tPow().field(), df.field());

    return tPow;
}


template<class Type, class GeoMesh, int r>
tmp<DimensionedField<typename powProduct<Type, r>::type, GeoMesh> >
pow
(
    const tmp<DimensionedField<Type, GeoMesh> >& tdf,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;

    const DimensionedField<Type, GeoMesh>& df = tdf();

    tmp<DimensionedField<powProductType, GeoMesh> > tPow =
        reuseTmpDimensionedField<powProductType, Type, GeoMesh>::New
        (
            tdf,
            "pow(" + df.name() + ',' + name(r) + ')',
            pow(df.dimensions(), r)
        );

    pow<Type, r, GeoMesh>(tPow().field(), df.field());

    reuseTmpDimensionedField<powProductType, Type, GeoMesh>::clear(tdf);

    return tPow;
}

template<class Type, class GeoMesh>
tmp<DimensionedField<typename outerProduct<Type, Type>::type, GeoMesh> >
sqr(const DimensionedField<Type, GeoMesh>& df)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;

    tmp<DimensionedField<outerProductType, GeoMesh> > tSqr
    (
        new DimensionedField<outerProductType, GeoMesh>
        (
            IOobject
            (
                "sqr(" + df.name() + ')',
                df.instance(),
                df.db()
            ),
            df.mesh(),
            sqr(df.dimensions())
        )
    );

    sqr(tSqr().field(), df.field());

    return tSqr;
}

template<class Type, class GeoMesh>
tmp<DimensionedField<typename outerProduct<Type, Type>::type, GeoMesh> >
sqr(const tmp<DimensionedField<Type, GeoMesh> >& tdf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;

    const DimensionedField<Type, GeoMesh>& df = tdf();

    tmp<DimensionedField<outerProductType, GeoMesh> > tSqr =
        reuseTmpDimensionedField<outerProductType, Type, GeoMesh>::New
        (
            tdf,
            "sqr(" + df.name() + ')',
            sqr(df.dimensions())
        );

    sqr(tSqr().field(), df.field());

    reuseTmpDimensionedField<outerProductType, Type, GeoMesh>::clear(tdf);

    return tSqr;
}


template<class Type, class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > magSqr
(
    const DimensionedField<Type, GeoMesh>& df
)
{
    tmp<DimensionedField<scalar, GeoMesh> > tMagSqr
    (
        new DimensionedField<scalar, GeoMesh>
        (
            IOobject
            (
                "magSqr(" + df.name() + ')',
                df.instance(),
                df.db()
            ),
            df.mesh(),
            sqr(df.dimensions())
        )
    );

    magSqr(tMagSqr().field(), df.field());

    return tMagSqr;
}

template<class Type, class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > magSqr
(
    const tmp<DimensionedField<Type, GeoMesh> >& tdf
)
{
    const DimensionedField<Type, GeoMesh>& df = tdf();

    tmp<DimensionedField<scalar, GeoMesh> > tMagSqr =
        reuseTmpDimensionedField<scalar, Type, GeoMesh>::New
        (
            tdf,
            "magSqr(" + df.name() + ')',
            sqr(df.dimensions())
        );

    magSqr(tMagSqr().field(), df.field());

    reuseTmpDimensionedField<scalar, Type, GeoMesh>::clear(tdf);

    return tMagSqr;
}


template<class Type, class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > mag
(
    const DimensionedField<Type, GeoMesh>& df
)
{
    tmp<DimensionedField<scalar, GeoMesh> > tMag
    (
        new DimensionedField<scalar, GeoMesh>
        (
            IOobject
            (
                "mag(" + df.name() + ')',
                df.instance(),
                df.db()
            ),
            df.mesh(),
            df.dimensions()
        )
    );

    mag(tMag().field(), df.field());

    return tMag;
}

template<class Type, class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > mag
(
    const tmp<DimensionedField<Type, GeoMesh> >& tdf
)
{
    const DimensionedField<Type, GeoMesh>& df = tdf();

    tmp<DimensionedField<scalar, GeoMesh> > tMag =
        reuseTmpDimensionedField<scalar, Type, GeoMesh>::New
        (
            tdf,
            "mag(" + df.name() + ')',
            df.dimensions()
        );

    mag(tMag().field(), df.field());

    reuseTmpDimensionedField<scalar, Type, GeoMesh>::clear(tdf);

    return tMag;
}


template<class Type, class GeoMesh>
tmp
<
    DimensionedField
        <typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh>
>
cmptAv(const DimensionedField<Type, GeoMesh>& df)
{
    typedef typename DimensionedField<Type, GeoMesh>::cmptType cmptType;

    tmp<DimensionedField<cmptType, GeoMesh> > CmptAv
    (
        new DimensionedField<scalar, GeoMesh>
        (
            IOobject
            (
                "cmptAv(" + df.name() + ')',
                df.instance(),
                df.db()
            ),
            df.mesh(),
            df.dimensions()
        )
    );

    cmptAv(CmptAv().field(), df.field());

    return CmptAv;
}

template<class Type, class GeoMesh>
tmp
<
    DimensionedField
        <typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh>
>
cmptAv(const tmp<DimensionedField<Type, GeoMesh> >& tdf)
{
    typedef typename DimensionedField<Type, GeoMesh>::cmptType
        cmptType;

    const DimensionedField<Type, GeoMesh>& df = tdf();

    tmp<DimensionedField<cmptType, GeoMesh> > CmptAv =
        reuseTmpDimensionedField<cmptType, Type, GeoMesh>::New
        (
            tdf,
            "cmptAv(" + df.name() + ')',
            df.dimensions()
        );

    cmptAv(CmptAv().field(), df.field());

    reuseTmpDimensionedField<cmptType, Type, GeoMesh>::clear(tdf);

    return CmptAv;
}

#define UNARY_REDUCTION_FUNCTION(returnType, func, dfunc)                     \
                                                                              \
template<class Type, class GeoMesh>                                           \
dimensioned<returnType> func                                                  \
(                                                                             \
    const DimensionedField<Type, GeoMesh>& df                                 \
)                                                                             \
{                                                                             \
    return dimensioned<Type>                                                  \
    (                                                                         \
        #func "(" + df.name() + ')',                                          \
        df.dimensions(),                                                      \
        dfunc(df.field())                                                     \
    );                                                                        \
}                                                                             \
                                                                              \
template<class Type, class GeoMesh>                                           \
dimensioned<returnType> func                                                  \
(                                                                             \
    const tmp<DimensionedField<Type, GeoMesh> >& tdf1                         \
)                                                                             \
{                                                                             \
    dimensioned<returnType> res = func(tdf1());                               \
    tdf1.clear();                                                             \
    return res;                                                               \
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

#define PRODUCT_OPERATOR(product, op, opFunc)                                 \
                                                                              \
template<class Type1, class Type2, class GeoMesh>                             \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh> >         \
operator op                                                                   \
(                                                                             \
    const DimensionedField<Type1, GeoMesh>& df1,                              \
    const DimensionedField<Type2, GeoMesh>& df2                               \
)                                                                             \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<DimensionedField<productType, GeoMesh> > tRes                         \
    (                                                                         \
        new DimensionedField<productType, GeoMesh>                            \
        (                                                                     \
            IOobject                                                          \
            (                                                                 \
                '(' + df1.name() + #op + df2.name() + ')',                    \
                df1.instance(),                                               \
                df1.db()                                                      \
            ),                                                                \
            df1.mesh(),                                                       \
            df1.dimensions() op df2.dimensions()                              \
        )                                                                     \
    );                                                                        \
                                                                              \
    Foam::opFunc(tRes().field(), df1.field(), df2.field());                   \
                                                                              \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Type1, class Type2, class GeoMesh>                             \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh> >         \
operator op                                                                   \
(                                                                             \
    const DimensionedField<Type1, GeoMesh>& df1,                              \
    const tmp<DimensionedField<Type2, GeoMesh> >& tdf2                        \
)                                                                             \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
                                                                              \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                     \
                                                                              \
    tmp<DimensionedField<productType, GeoMesh> > tRes =                       \
        reuseTmpDimensionedField<productType, Type2, GeoMesh>::New            \
        (                                                                     \
            tdf2,                                                             \
            '(' + df1.name() + #op + df2.name() + ')',                        \
            df1.dimensions() op df2.dimensions()                              \
        );                                                                    \
                                                                              \
    Foam::opFunc(tRes().field(), df1.field(), df2.field());                   \
                                                                              \
    reuseTmpDimensionedField<productType, Type2, GeoMesh>::clear(tdf2);       \
                                                                              \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Type1, class Type2, class GeoMesh>                             \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh> >         \
operator op                                                                   \
(                                                                             \
    const tmp<DimensionedField<Type1, GeoMesh> >& tdf1,                       \
    const DimensionedField<Type2, GeoMesh>& df2                               \
)                                                                             \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                     \
                                                                              \
    tmp<DimensionedField<productType, GeoMesh> > tRes =                       \
        reuseTmpDimensionedField<productType, Type1, GeoMesh>::New            \
        (                                                                     \
            tdf1,                                                             \
            '(' + df1.name() + #op + df2.name() + ')',                        \
            df1.dimensions() op df2.dimensions()                              \
        );                                                                    \
                                                                              \
    Foam::opFunc(tRes().field(), df1.field(), df2.field());                   \
                                                                              \
    reuseTmpDimensionedField<productType, Type1, GeoMesh>::clear(tdf1);       \
                                                                              \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Type1, class Type2, class GeoMesh>                             \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh> >         \
operator op                                                                   \
(                                                                             \
    const tmp<DimensionedField<Type1, GeoMesh> >& tdf1,                       \
    const tmp<DimensionedField<Type2, GeoMesh> >& tdf2                        \
)                                                                             \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                     \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                     \
                                                                              \
    tmp<DimensionedField<productType, GeoMesh> > tRes =                       \
        reuseTmpTmpDimensionedField                                           \
        <productType, Type1, Type1, Type2, GeoMesh>::New                      \
        (                                                                     \
            tdf1,                                                             \
            tdf2,                                                             \
            '(' + df1.name() + #op + df2.name() + ')',                        \
            df1.dimensions() op df2.dimensions()                              \
        );                                                                    \
                                                                              \
    Foam::opFunc(tRes().field(), df1.field(), df2.field());                   \
                                                                              \
    reuseTmpTmpDimensionedField                                               \
        <productType, Type1, Type1, Type2, GeoMesh>::clear(tdf1, tdf2);       \
                                                                              \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Form, class Type, class GeoMesh>                               \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const DimensionedField<Type, GeoMesh>& df1,                               \
    const dimensioned<Form>& dvs                                              \
)                                                                             \
{                                                                             \
    typedef typename product<Type, Form>::type productType;                   \
                                                                              \
    tmp<DimensionedField<productType, GeoMesh> > tRes                         \
    (                                                                         \
        new DimensionedField<productType, GeoMesh>                            \
        (                                                                     \
            IOobject                                                          \
            (                                                                 \
                '(' + df1.name() + #op + dvs.name() + ')',                    \
                df1.instance(),                                               \
                df1.db()                                                      \
            ),                                                                \
            df1.mesh(),                                                       \
            df1.dimensions() op dvs.dimensions()                              \
        )                                                                     \
    );                                                                        \
                                                                              \
    Foam::opFunc(tRes().field(), df1.field(), dvs.value());                   \
                                                                              \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type, class GeoMesh>        \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const DimensionedField<Type, GeoMesh>& df1,                               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
)                                                                             \
{                                                                             \
    return df1 op dimensioned<Form>(static_cast<const Form&>(vs));            \
}                                                                             \
                                                                              \
                                                                              \
template<class Form, class Type, class GeoMesh>                               \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const tmp<DimensionedField<Type, GeoMesh> >& tdf1,                        \
    const dimensioned<Form>& dvs                                              \
)                                                                             \
{                                                                             \
    typedef typename product<Type, Form>::type productType;                   \
                                                                              \
    const DimensionedField<Type, GeoMesh>& df1 = tdf1();                      \
                                                                              \
    tmp<DimensionedField<productType, GeoMesh> > tRes =                       \
        reuseTmpDimensionedField<productType, Type, GeoMesh>::New             \
        (                                                                     \
            tdf1,                                                             \
            '(' + df1.name() + #op + dvs.name() + ')',                        \
            df1.dimensions() op dvs.dimensions()                              \
        );                                                                    \
                                                                              \
    Foam::opFunc(tRes().field(), df1.field(), dvs.value());                   \
                                                                              \
    reuseTmpDimensionedField<productType, Type, GeoMesh>::clear(tdf1);        \
                                                                              \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type, class GeoMesh>        \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const tmp<DimensionedField<Type, GeoMesh> >& tdf1,                        \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
)                                                                             \
{                                                                             \
    return tdf1 op dimensioned<Form>(static_cast<const Form&>(vs));           \
}                                                                             \
                                                                              \
                                                                              \
template<class Form, class Type, class GeoMesh>                               \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const dimensioned<Form>& dvs,                                             \
    const DimensionedField<Type, GeoMesh>& df1                                \
)                                                                             \
{                                                                             \
    typedef typename product<Form, Type>::type productType;                   \
    tmp<DimensionedField<productType, GeoMesh> > tRes                         \
    (                                                                         \
        new DimensionedField<productType, GeoMesh>                            \
        (                                                                     \
            IOobject                                                          \
            (                                                                 \
                '(' + dvs.name() + #op + df1.name() + ')',                    \
                df1.instance(),                                               \
                df1.db()                                                      \
            ),                                                                \
            df1.mesh(),                                                       \
            dvs.dimensions() op df1.dimensions()                              \
        )                                                                     \
    );                                                                        \
                                                                              \
    Foam::opFunc(tRes().field(), dvs.value(), df1.field());                   \
                                                                              \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type, class GeoMesh>        \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const DimensionedField<Type, GeoMesh>& df1                                \
)                                                                             \
{                                                                             \
    return dimensioned<Form>(static_cast<const Form&>(vs)) op df1;            \
}                                                                             \
                                                                              \
template<class Form, class Type, class GeoMesh>                               \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const dimensioned<Form>& dvs,                                             \
    const tmp<DimensionedField<Type, GeoMesh> >& tdf1                         \
)                                                                             \
{                                                                             \
    typedef typename product<Form, Type>::type productType;                   \
                                                                              \
    const DimensionedField<Type, GeoMesh>& df1 = tdf1();                      \
                                                                              \
    tmp<DimensionedField<productType, GeoMesh> > tRes =                       \
        reuseTmpDimensionedField<productType, Type, GeoMesh>::New             \
        (                                                                     \
            tdf1,                                                             \
            '(' + dvs.name() + #op + df1.name() + ')',                        \
            dvs.dimensions() op df1.dimensions()                              \
        );                                                                    \
                                                                              \
    Foam::opFunc(tRes().field(), dvs.value(), df1.field());                   \
                                                                              \
    reuseTmpDimensionedField<productType, Type, GeoMesh>::clear(tdf1);        \
                                                                              \
    return tRes;                                                              \
}                                                                             \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type, class GeoMesh>        \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const tmp<DimensionedField<Type, GeoMesh> >& tdf1                         \
)                                                                             \
{                                                                             \
    return dimensioned<Form>(static_cast<const Form&>(vs)) op tdf1;           \
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
