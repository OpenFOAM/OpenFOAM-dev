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

#include "FieldM.H"
#include "FieldFieldReuseFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type, Func)                                \
                                                                              \
TEMPLATE                                                                      \
void Func                                                                     \
(                                                                             \
    FieldField<Field, ReturnType>& res,                                       \
    const FieldField<Field, Type>& f                                          \
)                                                                             \
{                                                                             \
    forAll(res, i)                                                            \
    {                                                                         \
        Func(res[i], f[i]);                                                   \
    }                                                                         \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const FieldField<Field, Type>& f                                          \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        FieldField<Field, ReturnType>::NewCalculatedType(f)                   \
    );                                                                        \
    Func(tRes(), f);                                                          \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf                                   \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpFieldField<Field, Type, Type>::New(tf)                        \
    );                                                                        \
    Func(tRes(), tf());                                                       \
    reuseTmpFieldField<Field, Type, Type>::clear(tf);                         \
    return tRes;                                                              \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                          \
                                                                              \
TEMPLATE                                                                      \
void OpFunc                                                                   \
(                                                                             \
    FieldField<Field, ReturnType>& res,                                       \
    const FieldField<Field, Type>& f                                          \
)                                                                             \
{                                                                             \
    forAll(res, i)                                                            \
    {                                                                         \
        OpFunc(res[i], f[i]);                                                 \
    }                                                                         \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const FieldField<Field, Type>& f                                          \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        FieldField<Field, Type>::NewCalculatedType(f)                         \
    );                                                                        \
    OpFunc(tRes(), f);                                                        \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf                                   \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpFieldField<Field, Type, Type>::New(tf)                        \
    );                                                                        \
    OpFunc(tRes(), tf());                                                     \
    reuseTmpFieldField<Field, Type, Type>::clear(tf);                         \
    return tRes;                                                              \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                       \
                                                                              \
TEMPLATE                                                                      \
void Func                                                                     \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        Func(f[i], f1[i], f2[i]);                                             \
    }                                                                         \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        FieldField<Field, Type1>::NewCalculatedType(f1)                       \
    );                                                                        \
    Func(tRes(), f1, f2);                                                     \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2)                \
    );                                                                        \
    Func(tRes(), f1, tf2());                                                  \
    reuseTmpFieldField<Field, ReturnType, Type2>::clear(tf2);                 \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1)                \
    );                                                                        \
    Func(tRes(), tf1(), f2);                                                  \
    reuseTmpFieldField<Field, ReturnType, Type1>::clear(tf1);                 \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::       \
            New(tf1, tf2)                                                     \
    );                                                                        \
    Func(tRes(), tf1(), tf2());                                               \
    reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::           \
        clear(tf1, tf2);                                                      \
    return tRes;                                                              \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)               \
                                                                              \
TEMPLATE                                                                      \
void Func                                                                     \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const FieldField<Field, Type1>& f1,                                       \
    const Type2& s                                                            \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        Func(f[i], f1[i], s);                                                 \
    }                                                                         \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const Type2& s                                                            \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        FieldField<Field, Type1>::NewCalculatedType(f1)                       \
    );                                                                        \
    Func(tRes(), f1, s);                                                      \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const Type2& s                                                            \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1)                \
    );                                                                        \
    Func(tRes(), tf1(), s);                                                   \
    reuseTmpFieldField<Field, ReturnType, Type1>::clear(tf1);                 \
    return tRes;                                                              \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)               \
                                                                              \
TEMPLATE                                                                      \
void Func                                                                     \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const Type1& s,                                                           \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        Func(f[i], s, f2[i]);                                                 \
    }                                                                         \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const Type1& s,                                                           \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        FieldField<Field, Type2>::NewCalculatedType(f2)                       \
    );                                                                        \
    Func(tRes(), s, f2);                                                      \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const Type1& s,                                                           \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2)                \
    );                                                                        \
    Func(tRes(), s, tf2());                                                   \
    reuseTmpFieldField<Field, ReturnType, Type2>::clear(tf2);                 \
    return tRes;                                                              \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                  \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                 \
                                                                              \
TEMPLATE                                                                      \
void OpFunc                                                                   \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        OpFunc(f[i], f1[i], f2[i]);                                           \
    }                                                                         \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        FieldField<Field, ReturnType>::NewCalculatedType(f1)                  \
    );                                                                        \
    OpFunc(tRes(), f1, f2);                                                   \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2)                \
    );                                                                        \
    OpFunc(tRes(), f1, tf2());                                                \
    reuseTmpFieldField<Field, ReturnType, Type2>::clear(tf2);                 \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1)                \
    );                                                                        \
    OpFunc(tRes(), tf1(), f2);                                                \
    reuseTmpFieldField<Field, ReturnType, Type1>::clear(tf1);                 \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::       \
            New(tf1, tf2)                                                     \
    );                                                                        \
    OpFunc(tRes(), tf1(), tf2());                                             \
    reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::           \
        clear(tf1, tf2);                                                      \
    return tRes;                                                              \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)         \
                                                                              \
TEMPLATE                                                                      \
void OpFunc                                                                   \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const Type1& s,                                                           \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        OpFunc(f[i], s, f2[i]);                                               \
    }                                                                         \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const Type1& s,                                                           \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        FieldField<Field, Type2>::NewCalculatedType(f2)                       \
    );                                                                        \
    OpFunc(tRes(), s, f2);                                                    \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const Type1& s,                                                           \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2)                \
    );                                                                        \
    OpFunc(tRes(), s, tf2());                                                 \
    reuseTmpFieldField<Field, ReturnType, Type2>::clear(tf2);                 \
    return tRes;                                                              \
}


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)         \
                                                                              \
TEMPLATE                                                                      \
void OpFunc                                                                   \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const FieldField<Field, Type1>& f1,                                       \
    const Type2& s                                                            \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        OpFunc(f[i], f1[i], s);                                               \
    }                                                                         \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const Type2& s                                                            \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        FieldField<Field, Type1>::NewCalculatedType(f1)                       \
    );                                                                        \
    OpFunc(tRes(), f1, s);                                                    \
    return tRes;                                                              \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const Type2& s                                                            \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, ReturnType> > tRes                                  \
    (                                                                         \
        reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1)                \
    );                                                                        \
    OpFunc(tRes(), tf1(), s);                                                 \
    reuseTmpFieldField<Field, ReturnType, Type1>::clear(tf1);                 \
    return tRes;                                                              \
}


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)            \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)

// ************************************************************************* //
