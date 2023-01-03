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
#include "GeometricFieldFunctionsM.H"

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
);

TEMPLATE2
void T
(
     GeometricField<Type, GeoMesh, PrimitiveField1>& gf,
     const GeometricField<Type, GeoMesh, PrimitiveField2>& gf1
);

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
);

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
);

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
);

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
);

TEMPLATE
tmp<GeometricField<typename outerProduct<Type, Type>::type, GeoMesh, Field>>
sqr(const GeometricField<Type, GeoMesh, PrimitiveField>& gf);

TEMPLATE
tmp<GeometricField<typename outerProduct<Type, Type>::type, GeoMesh, Field>>
sqr(const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf);

TEMPLATE2
void magSqr
(
    GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf
);

TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> magSqr
(
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf
);

TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> magSqr
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
);

TEMPLATE2
void mag
(
    GeometricField<scalar, GeoMesh, PrimitiveField1>& gsf,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf
);

TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> mag
(
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf
);

TEMPLATE
tmp<GeometricField<scalar, GeoMesh, Field>> mag
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
);

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
);

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
cmptAv(const GeometricField<Type, GeoMesh, PrimitiveField>& gf);

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
cmptAv(const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf);

UNARY_FUNCTION(Type, Type, cmptMag, cmptMag);


#define UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(returnType, func, gFunc)        \
                                                                               \
TEMPLATE                                                                       \
dimensioned<returnType> func                                                   \
(                                                                              \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf                    \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
dimensioned<returnType> func                                                   \
(                                                                              \
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf1             \
);

UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(Type, max, gMax)
UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(Type, min, gMin)

#undef UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY


#define UNARY_REDUCTION_FUNCTION(returnType, func, gFunc)                      \
                                                                               \
TEMPLATE                                                                       \
dimensioned<returnType> func                                                   \
(                                                                              \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf                    \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
dimensioned<returnType> func                                                   \
(                                                                              \
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf1             \
);

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
tmp<GeometricField<typename product<Type1, Type2>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField1>& gf1,                \
    const GeometricField<Type2, GeoMesh, PrimitiveField2>& gf2                 \
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
tmp<GeometricField<typename product<Type1, Type2>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField1>& gf1,                \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField2>>& tgf2           \
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
tmp<GeometricField<typename product<Type1, Type2>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField1>>& tgf1,          \
    const GeometricField<Type2, GeoMesh, PrimitiveField2>& gf2                 \
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
tmp<GeometricField<typename product<Type1, Type2>::type, GeoMesh, Field>>      \
operator op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField1>>& tgf1,          \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField2>>& tgf2           \
);                                                                             \
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
);                                                                             \
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
tmp<GeometricField<typename product<Form, Type>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf1,                  \
    const VectorSpace<Form, Cmpt, nCmpt>& vs                                   \
);                                                                             \
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
tmp<GeometricField<typename product<Form, Type>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf1,            \
    const VectorSpace<Form, Cmpt, nCmpt>& vs                                   \
);                                                                             \
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
);                                                                             \
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
tmp<GeometricField<typename product<Form, Type>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const VectorSpace<Form, Cmpt, nCmpt>& vs,                                  \
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf1                   \
);                                                                             \
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
tmp<GeometricField<typename product<Form, Type>::type, GeoMesh, Field>>        \
operator op                                                                    \
(                                                                              \
    const VectorSpace<Form, Cmpt, nCmpt>& vs,                                  \
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf1             \
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
