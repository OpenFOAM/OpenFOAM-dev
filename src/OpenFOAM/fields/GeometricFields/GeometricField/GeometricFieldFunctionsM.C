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

#include "GeometricFieldReuseFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type1, Func, Dfunc)                         \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& res,                      \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1                      \
)                                                                              \
{                                                                              \
    Foam::Func(res.primitiveFieldRef(), gf1.primitiveField());                 \
    Foam::Func(res.boundaryFieldRef(), gf1.boundaryField());                   \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1                      \
)                                                                              \
{                                                                              \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        GeometricField<ReturnType, PatchField, GeoMesh>::New                   \
        (                                                                      \
            #Func "(" + gf1.name() + ')',                                      \
            gf1.mesh(),                                                        \
            Dfunc(gf1.dimensions())                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1);                                               \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1                \
)                                                                              \
{                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf1,                                                              \
            #Func "(" + gf1.name() + ')',                                      \
            Dfunc(gf1.dimensions())                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1);                                               \
                                                                               \
    tgf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc, Dfunc)                   \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& res,                      \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1                      \
)                                                                              \
{                                                                              \
    Foam::OpFunc(res.primitiveFieldRef(), gf1.primitiveField());               \
    Foam::OpFunc(res.boundaryFieldRef(), gf1.boundaryField());                 \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1                      \
)                                                                              \
{                                                                              \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        GeometricField<ReturnType, PatchField, GeoMesh>::New                   \
        (                                                                      \
            #Op + gf1.name(),                                                  \
            gf1.mesh(),                                                        \
            Dfunc(gf1.dimensions())                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref(), gf1);                                             \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1                \
)                                                                              \
{                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf1,                                                              \
            #Op + gf1.name(),                                                  \
            Dfunc(gf1.dimensions())                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref(), gf1);                                             \
                                                                               \
    tgf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& res,                      \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    Foam::Func                                                                 \
    (                                                                          \
        res.primitiveFieldRef(),                                               \
        gf1.primitiveField(),                                                  \
        gf2.primitiveField()                                                   \
    );                                                                         \
    Foam::Func                                                                 \
    (                                                                          \
        res.boundaryFieldRef(),                                                \
        gf1.boundaryField(),                                                   \
        gf2.boundaryField()                                                    \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        GeometricField<ReturnType, PatchField, GeoMesh>::New                   \
        (                                                                      \
            #Func "(" + gf1.name() + ',' + gf2.name() + ')',                   \
            gf1.mesh(),                                                        \
            Func(gf1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, gf2);                                          \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2 = tgf2();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf2,                                                              \
            #Func "(" + gf1.name() + ',' + gf2.name() + ')',                   \
            Func(gf1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, gf2);                                          \
                                                                               \
    tgf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf1,                                                              \
            #Func "(" + gf1.name() + ',' + gf2.name() + ')',                   \
            Func(gf1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, gf2);                                          \
                                                                               \
    tgf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();            \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2 = tgf2();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpTmpGeometricField                                              \
        <ReturnType, Type1, Type2, PatchField, GeoMesh>::New                   \
        (                                                                      \
            tgf1,                                                              \
            tgf2,                                                              \
            #Func "(" + gf1.name() + ',' + gf2.name() + ')',                   \
            Func(gf1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, gf2);                                          \
                                                                               \
    tgf1.clear();                                                              \
    tgf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& res,                      \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    Foam::Func(res.primitiveFieldRef(), dt1.value(), gf2.primitiveField());    \
    Foam::Func(res.boundaryFieldRef(), dt1.value(), gf2.boundaryField());      \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        GeometricField<ReturnType, PatchField, GeoMesh>::New                   \
        (                                                                      \
            #Func "(" + dt1.name() + ',' + gf2.name() + ')',                   \
            gf2.mesh(),                                                        \
            Func(dt1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), dt1, gf2);                                          \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const Type1& t1,                                                           \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    return Func(dimensioned<Type1>(t1), gf2);                                  \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2 = tgf2();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf2,                                                              \
            #Func "(" + dt1.name() + gf2.name() + ',' + ')',                   \
            Func(dt1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), dt1, gf2);                                          \
                                                                               \
    tgf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const Type1& t1,                                                           \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    return Func(dimensioned<Type1>(t1), tgf2);                                 \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& res,                      \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    Foam::Func(res.primitiveFieldRef(), gf1.primitiveField(), dt2.value());    \
    Foam::Func(res.boundaryFieldRef(), gf1.boundaryField(), dt2.value());      \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        GeometricField<ReturnType, PatchField, GeoMesh>::New                   \
        (                                                                      \
            #Func "(" + gf1.name() + ',' + dt2.name() + ')',                   \
            gf1.mesh(),                                                        \
            Func(gf1.dimensions(), dt2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, dt2);                                          \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return Func(gf1, dimensioned<Type2>(t2));                                  \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf1,                                                              \
            #Func "(" + gf1.name() + ',' + dt2.name() + ')',                   \
            Func(gf1.dimensions(), dt2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, dt2);                                          \
                                                                               \
    tgf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return Func(tgf1, dimensioned<Type2>(t2));                                 \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpName, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& res,                      \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    Foam::OpFunc                                                               \
    (res.primitiveFieldRef(), gf1.primitiveField(), gf2.primitiveField());     \
    Foam::OpFunc                                                               \
    (res.boundaryFieldRef(), gf1.boundaryField(), gf2.boundaryField());        \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        GeometricField<ReturnType, PatchField, GeoMesh>::New                   \
        (                                                                      \
            '(' + gf1.name() + OpName + gf2.name() + ')',                      \
            gf1.mesh(),                                                        \
            gf1.dimensions() Op gf2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref(), gf1, gf2);                                        \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2 = tgf2();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf2,                                                              \
            '(' + gf1.name() + OpName + gf2.name() + ')',                      \
            gf1.dimensions() Op gf2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref(), gf1, gf2);                                        \
                                                                               \
    tgf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf1,                                                              \
            '(' + gf1.name() + OpName + gf2.name() + ')',                      \
            gf1.dimensions() Op gf2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref(), gf1, gf2);                                        \
                                                                               \
    tgf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();            \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2 = tgf2();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpTmpGeometricField                                              \
        <ReturnType, Type1, Type2, PatchField, GeoMesh>::New                   \
        (                                                                      \
            tgf1,                                                              \
            tgf2,                                                              \
            '(' + gf1.name() + OpName + gf2.name() + ')',                      \
            gf1.dimensions() Op gf2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref(), gf1, gf2);                                        \
                                                                               \
    tgf1.clear();                                                              \
    tgf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc)  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& res,                      \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    Foam::OpFunc(res.primitiveFieldRef(), dt1.value(), gf2.primitiveField());  \
    Foam::OpFunc(res.boundaryFieldRef(), dt1.value(), gf2.boundaryField());    \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        GeometricField<ReturnType, PatchField, GeoMesh>::New                   \
        (                                                                      \
            '(' + dt1.name() + OpName + gf2.name() + ')',                      \
            gf2.mesh(),                                                        \
            dt1.dimensions() Op gf2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref(), dt1, gf2);                                        \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const Type1& t1,                                                           \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    return dimensioned<Type1>(t1) Op gf2;                                      \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2 = tgf2();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf2,                                                              \
            '(' + dt1.name() + OpName + gf2.name() + ')',                      \
            dt1.dimensions() Op gf2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref(), dt1, gf2);                                        \
                                                                               \
    tgf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const Type1& t1,                                                           \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    return dimensioned<Type1>(t1) Op tgf2;                                     \
}


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc)  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& res,                      \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    Foam::OpFunc(res.primitiveFieldRef(), gf1.primitiveField(), dt2.value());  \
    Foam::OpFunc(res.boundaryFieldRef(), gf1.boundaryField(), dt2.value());    \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        GeometricField<ReturnType, PatchField, GeoMesh>::New                   \
        (                                                                      \
            '(' + gf1.name() + OpName + dt2.name() + ')',                      \
            gf1.mesh(),                                                        \
            gf1.dimensions() Op dt2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref(), gf1, dt2);                                        \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return gf1 Op dimensioned<Type2>(t2);                                      \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf1,                                                              \
            '(' + gf1.name() + OpName + dt2.name() + ')',                      \
            gf1.dimensions() Op dt2.dimensions()                               \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::OpFunc(tRes.ref(), gf1, dt2);                                        \
                                                                               \
    tgf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return tgf1 Op dimensioned<Type2>(t2);                                     \
}


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpName, OpFunc)     \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc)      \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
