/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
#include "FieldReuseFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type, Func)                                 \
                                                                               \
TEMPLATE                                                                       \
void Func(Field<ReturnType>& res, const UList<Type>& f)                        \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_F(ReturnType, res, =, ::Foam::Func, Type, f)            \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                           \
                                                                               \
TEMPLATE                                                                       \
void OpFunc(Field<ReturnType>& res, const UList<Type>& f)                      \
{                                                                              \
    TFOR_ALL_F_OP_OP_F(ReturnType, res, =, Op, Type, f)                        \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_F_F                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, f1, Type2, f2                 \
    )                                                                          \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const Type1& s1,                                                           \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_S_F                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, s1, Type2, f2                 \
    )                                                                          \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const UList<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_F_S                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, f1, Type2, s2                 \
    )                                                                          \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_F_OP_F(ReturnType, res, =, Type1, f1, Op, Type2, f2)         \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const Type1& s1,                                                           \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_S_OP_F(ReturnType, res, =, Type1, s1, Op, Type2, f2)         \
}


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const UList<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_F_OP_S(ReturnType, res, =, Type1, f1, Op, Type2, s2)         \
}


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)


// ************************************************************************* //
