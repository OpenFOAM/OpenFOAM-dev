/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
#include "DimensionedFieldFunctionsM.H"

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
);

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
);

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
sqr(const DimensionedField<Type, GeoMesh, PrimitiveField>& df);

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
sqr(const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf);

template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> magSqr
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df
);

template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> magSqr
(
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf
);

template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> mag
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df
);

template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<scalar, GeoMesh, Field>> mag
(
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf
);

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
cmptAv(const DimensionedField<Type, GeoMesh, PrimitiveField>& df);

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
cmptAv(const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf);

UNARY_FUNCTION(Type, Type, cmptMag, cmptMag);


#define UNARY_REDUCTION_FUNCTION(returnType, func, dfunc)                      \
                                                                               \
template<class Type, class GeoMesh, template<class> class PrimitiveField>      \
dimensioned<returnType> func                                                   \
(                                                                              \
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df                  \
);                                                                             \
template<class Type, class GeoMesh, template<class> class PrimitiveField>      \
dimensioned<returnType> func                                                   \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf1           \
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
);                                                                             \
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
);                                                                             \
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
);                                                                             \
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
);                                                                             \
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
);                                                                             \
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
);                                                                             \
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
);                                                                             \
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
);                                                                             \
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
);                                                                             \
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
);                                                                             \
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
);                                                                             \
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
