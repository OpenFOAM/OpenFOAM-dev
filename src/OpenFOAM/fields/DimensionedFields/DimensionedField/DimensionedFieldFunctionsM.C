/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type1, Func, Dfunc)                         \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1                                \
)                                                                              \
{                                                                              \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        DimensionedField<ReturnType, GeoMesh>::New                             \
        (                                                                      \
            #Func "(" + df1.name() + ')',                                      \
            df1.mesh(),                                                        \
            Dfunc(df1.dimensions())                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field());                                     \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1                          \
)                                                                              \
{                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type1, GeoMesh>::New              \
        (                                                                      \
            tdf1,                                                              \
            #Func "(" + df1.name() + ')',                                      \
            Dfunc(df1.dimensions())                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field());                                     \
                                                                               \
    tdf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc, Dfunc)                   \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1                                \
)                                                                              \
{                                                                              \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        DimensionedField<ReturnType, GeoMesh>::New                             \
        (                                                                      \
            #Op + df1.name(),                                                  \
            df1.mesh(),                                                        \
            Dfunc(df1.dimensions())                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref().field(), df1.field());                             \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1                          \
)                                                                              \
{                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type1, GeoMesh>::New              \
        (                                                                      \
            tdf1,                                                              \
            #Op + df1.name(),                                                  \
            Dfunc(df1.dimensions())                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref().field(), df1.field());                             \
                                                                               \
    tdf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        DimensionedField<ReturnType, GeoMesh>::New                             \
        (                                                                      \
            #Func "(" + df1.name() + ',' + df2.name() + ')',                   \
            df1.mesh(),                                                        \
            Func(df1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), df2.field());                        \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type2, GeoMesh>::New              \
        (                                                                      \
            tdf2,                                                              \
            #Func "(" + df1.name() + ',' + df2.name() + ')',                   \
            Func(df1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), df2.field());                        \
                                                                               \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type1, GeoMesh>::New              \
        (                                                                      \
            tdf1,                                                              \
            #Func "(" + df1.name() + ',' + df2.name() + ')',                   \
            Func(df1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), df2.field());                        \
                                                                               \
    tdf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpTmpDimensionedField<ReturnType, Type1, Type2, GeoMesh>::New    \
        (                                                                      \
            tdf1,                                                              \
            tdf2,                                                              \
            #Func "(" + df1.name() + ',' + df2.name() + ')',                   \
            Func(df1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), df2.field());                        \
                                                                               \
    tdf1.clear();                                                              \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        DimensionedField<ReturnType, GeoMesh>::New                             \
        (                                                                      \
            #Func "(" + dt1.name() + ',' + df2.name() + ')',                   \
            df2.mesh(),                                                        \
            Func(dt1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), dt1.value(), df2.field());                        \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const Type1& t1,                                                           \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    return Func(dimensioned<Type1>(t1), df2);                                  \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type2, GeoMesh>::New              \
        (                                                                      \
            tdf2,                                                              \
            #Func "(" + dt1.name() + ',' + df2.name() + ')',                   \
            Func(dt1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), dt1.value(), df2.field());                        \
                                                                               \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const Type1& t1,                                                           \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    return Func(dimensioned<Type2>(t1), tdf2);                                 \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        DimensionedField<ReturnType, GeoMesh>::New                             \
        (                                                                      \
            #Func "(" + df1.name() + ',' + dt2.name() + ')',                   \
            df1.mesh(),                                                        \
            Func(df1.dimensions(), dt2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), dt2.value());                        \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return Func(df1, dimensioned<Type2>(t2));                                  \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type1, GeoMesh>::New              \
        (                                                                      \
            tdf1,                                                              \
            #Func "(" + df1.name() + ',' + dt2.name() + ')',                   \
            Func(df1.dimensions(), dt2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), dt2.value());                        \
                                                                               \
    tdf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return Func(tdf1, dimensioned<Type2>(t2));                                 \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpName, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        DimensionedField<ReturnType, GeoMesh>::New                             \
        (                                                                      \
            '(' + df1.name() + OpName + df2.name() + ')',                      \
            df1.mesh(),                                                        \
            df1.dimensions() Op df2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref().field(), df1.field(), df2.field());                \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type2, GeoMesh>::New              \
        (                                                                      \
            tdf2,                                                              \
            '(' + df1.name() + OpName + df2.name() + ')',                      \
            df1.dimensions() Op df2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref().field(), df1.field(), df2.field());                \
                                                                               \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type1, GeoMesh>::New              \
        (                                                                      \
            tdf1,                                                              \
            '(' + df1.name() + OpName + df2.name() + ')',                      \
            df1.dimensions() Op df2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref().field(), df1.field(), df2.field());                \
                                                                               \
    tdf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpTmpDimensionedField<ReturnType, Type1, Type2, GeoMesh>::New    \
        (                                                                      \
            tdf1,                                                              \
            tdf2,                                                              \
            '(' + df1.name() + OpName + df2.name() + ')',                      \
            df1.dimensions() Op df2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref().field(), df1.field(), df2.field());                \
                                                                               \
    tdf1.clear();                                                              \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc) \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        DimensionedField<ReturnType, GeoMesh>::New                             \
        (                                                                      \
            '(' + dt1.name() + OpName + df2.name() + ')',                      \
            df2.mesh(),                                                        \
            dt1.dimensions() Op df2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref().field(), dt1.value(), df2.field());                \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const Type1& t1,                                                           \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    return dimensioned<Type1>(t1) Op df2;                                      \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type2, GeoMesh>::New              \
        (                                                                      \
            tdf2,                                                              \
            '(' + dt1.name() + OpName + df2.name() + ')',                      \
            dt1.dimensions() Op df2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref().field(), dt1.value(), tdf2().field());             \
                                                                               \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const Type1& t1,                                                           \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    return dimensioned<Type1>(t1) Op tdf2;                                     \
}


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc) \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        DimensionedField<ReturnType, GeoMesh>::New                             \
        (                                                                      \
            '(' + df1.name() + OpName + dt2.name() + ')',                      \
            df1.mesh(),                                                        \
            df1.dimensions() Op dt2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref().field(), df1.field(), dt2.value());                \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return df1 Op dimensioned<Type2>(t2);                                      \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type1, GeoMesh>::New              \
        (                                                                      \
            tdf1,                                                              \
            '(' + df1.name() + OpName + dt2.name() + ')',                      \
            df1.dimensions() Op dt2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref().field(), tdf1().field(), dt2.value());             \
                                                                               \
    tdf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> operator Op                         \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return tdf1 Op dimensioned<Type2>(t2);                                     \
}

#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpName, OpFunc)     \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc)      \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc)


// ************************************************************************* //
