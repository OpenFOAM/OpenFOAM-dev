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

#include "DimensionedScalarField.H"

#define TEMPLATE template<class Type, class GeoMesh>
#include "DimensionedFieldFunctionsM.H"

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
);

template<class Type, class GeoMesh, int r>
tmp<DimensionedField<typename powProduct<Type, r>::type, GeoMesh> >
pow
(
    const tmp<DimensionedField<Type, GeoMesh> >& tdf,
    typename powProduct<Type, r>::type
);

template<class Type, class GeoMesh>
tmp<DimensionedField<typename outerProduct<Type, Type>::type, GeoMesh> >
sqr(const DimensionedField<Type, GeoMesh>& df);

template<class Type, class GeoMesh>
tmp<DimensionedField<typename outerProduct<Type, Type>::type, GeoMesh> >
sqr(const tmp<DimensionedField<Type, GeoMesh> >& tdf);

template<class Type, class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > magSqr
(
    const DimensionedField<Type, GeoMesh>& df
);

template<class Type, class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > magSqr
(
    const tmp<DimensionedField<Type, GeoMesh> >& tdf
);

template<class Type, class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > mag
(
    const DimensionedField<Type, GeoMesh>& df
);

template<class Type, class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh> > mag
(
    const tmp<DimensionedField<Type, GeoMesh> >& tdf
);

template<class Type, class GeoMesh>
tmp
<
    DimensionedField
        <typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh>
>
cmptAv(const DimensionedField<Type, GeoMesh>& df);

template<class Type, class GeoMesh>
tmp
<
    DimensionedField
        <typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh>
>
cmptAv(const tmp<DimensionedField<Type, GeoMesh> >& tdf);


#define UNARY_REDUCTION_FUNCTION(returnType, func, dfunc)                     \
                                                                              \
template<class Type, class GeoMesh>                                           \
dimensioned<returnType> func                                                  \
(                                                                             \
    const DimensionedField<Type, GeoMesh>& df                                 \
);                                                                            \
template<class Type, class GeoMesh>                                           \
dimensioned<returnType> func                                                  \
(                                                                             \
    const tmp<DimensionedField<Type, GeoMesh> >& tdf1                         \
);

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
);                                                                            \
                                                                              \
template<class Type1, class Type2, class GeoMesh>                             \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh> >         \
operator op                                                                   \
(                                                                             \
    const DimensionedField<Type1, GeoMesh>& df1,                              \
    const tmp<DimensionedField<Type2, GeoMesh> >& tdf2                        \
);                                                                            \
                                                                              \
template<class Type1, class Type2, class GeoMesh>                             \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh> >         \
operator op                                                                   \
(                                                                             \
    const tmp<DimensionedField<Type1, GeoMesh> >& tdf1,                       \
    const DimensionedField<Type2, GeoMesh>& df2                               \
);                                                                            \
                                                                              \
template<class Type1, class Type2, class GeoMesh>                             \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh> >         \
operator op                                                                   \
(                                                                             \
    const tmp<DimensionedField<Type1, GeoMesh> >& tdf1,                       \
    const tmp<DimensionedField<Type2, GeoMesh> >& tdf2                        \
);                                                                            \
                                                                              \
template<class Form, class Type, class GeoMesh>                               \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const DimensionedField<Type, GeoMesh>& df1,                               \
    const dimensioned<Form>& dvs                                              \
);                                                                            \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type, class GeoMesh>        \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const DimensionedField<Type, GeoMesh>& df1,                               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
);                                                                            \
                                                                              \
template<class Form, class Type, class GeoMesh>                               \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const tmp<DimensionedField<Type, GeoMesh> >& tdf1,                        \
    const dimensioned<Form>& dvs                                              \
);                                                                            \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type, class GeoMesh>        \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const tmp<DimensionedField<Type, GeoMesh> >& tdf1,                        \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
);                                                                            \
                                                                              \
template<class Form, class Type, class GeoMesh>                               \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const dimensioned<Form>& dvs,                                             \
    const DimensionedField<Type, GeoMesh>& df1                                \
);                                                                            \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type, class GeoMesh>        \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const DimensionedField<Type, GeoMesh>& df1                                \
);                                                                            \
                                                                              \
template<class Form, class Type, class GeoMesh>                               \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const dimensioned<Form>& dvs,                                             \
    const tmp<DimensionedField<Type, GeoMesh> >& tdf1                         \
);                                                                            \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type, class GeoMesh>        \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh> >           \
operator op                                                                   \
(                                                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const tmp<DimensionedField<Type, GeoMesh> >& tdf1                         \
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
