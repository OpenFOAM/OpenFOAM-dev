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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type1, Func, Dfunc)                         \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField>& df                   \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField>>& tdf1            \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type1, Op, opFunc, Dfunc)                   \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField>& df1                  \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField>>& tdf1            \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE2                                                                      \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField1>& df1,                \
    const GeometricField<Type2, GeoMesh, PrimitiveField2>& df2                 \
);                                                                             \
                                                                               \
TEMPLATE2                                                                      \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField1>& df1,                \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField2>>& tdf2           \
);                                                                             \
                                                                               \
TEMPLATE2                                                                      \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField1>>& tdf1,          \
    const GeometricField<Type2, GeoMesh, PrimitiveField2>& df2                 \
);                                                                             \
                                                                               \
TEMPLATE2                                                                      \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField1>>& tdf1,          \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField2>>& tdf2           \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, GeoMesh, PrimitiveField>& df2                  \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const Type1& t1,                                                           \
    const GeometricField<Type2, GeoMesh, PrimitiveField>& df2                  \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField>>& tdf2            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const Type1& t1,                                                           \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField>>& tdf2            \
);


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField>& df1,                 \
    const dimensioned<Type2>& dt2                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField>& df1,                 \
    const Type2& t2                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField>>& tdf1,           \
    const dimensioned<Type2>& dt2                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> Func                           \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField>>& tdf2,           \
    const Type2& t2                                                            \
);


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpName, OpFunc)          \
                                                                               \
TEMPLATE2                                                                      \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField1>& df1,                \
    const GeometricField<Type2, GeoMesh, PrimitiveField2>& df2                 \
);                                                                             \
                                                                               \
TEMPLATE2                                                                      \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField1>& df1,                \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField2>>& tdf2           \
);                                                                             \
                                                                               \
TEMPLATE2                                                                      \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField1>>& tdf1,          \
    const GeometricField<Type2, GeoMesh, PrimitiveField2>& df2                 \
);                                                                             \
                                                                               \
TEMPLATE2                                                                      \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField1>>& tdf1,          \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField2>>& tdf2           \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc)  \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, GeoMesh, PrimitiveField>& df2                  \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const Type1& t1,                                                           \
    const GeometricField<Type2, GeoMesh, PrimitiveField>& df2                  \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField>>& tdf2            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const Type1& t1,                                                           \
    const tmp<GeometricField<Type2, GeoMesh, PrimitiveField>>& tdf2            \
);


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc)  \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField>& df1,                 \
    const dimensioned<Type2>& dt2                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const GeometricField<Type1, GeoMesh, PrimitiveField>& df1,                 \
    const Type2& t2                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField>>& tdf1,           \
    const dimensioned<Type2>& dt2                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, GeoMesh, Field>> operator Op                    \
(                                                                              \
    const tmp<GeometricField<Type1, GeoMesh, PrimitiveField>>& tdf1,           \
    const Type2& t2                                                            \
);


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpName, OpFunc)     \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc)      \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc)


// ************************************************************************* //
