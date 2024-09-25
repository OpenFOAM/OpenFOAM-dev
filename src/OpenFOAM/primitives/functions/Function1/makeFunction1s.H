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

Class
    Foam::Function1

Description
    Macros for creating standard Function1-s

\*---------------------------------------------------------------------------*/

#ifndef makeFunction1s_H
#define makeFunction1s_H

#include "None.H"
#include "Constant.H"
#include "Uniform.H"
#include "ZeroConstant.H"
#include "OneConstant.H"
#include "Polynomial1.H"
#include "Add.H"
#include "Repeat.H"
#include "Scale.H"
#include "Shift.H"
#include "Sine.H"
#include "Square.H"
#include "Table.H"
#include "UniformTable1.H"
#include "NonUniformTable1.H"
#include "CodedFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction1s(Type, nullArg)                                          \
                                                                               \
    template<>                                                                 \
    const char* const Foam::Tuple2<Foam::scalar, Type>::typeName               \
    (                                                                          \
        "Tuple2<scalar," #Type ">"                                             \
    );                                                                         \
                                                                               \
    defineFunction1(Type);                                                     \
                                                                               \
    namespace Function1s                                                       \
    {                                                                          \
        addFunction1(None, Type);                                              \
        addStreamConstructableFunction1(Constant, Type);                       \
        addStreamConstructableFunction1(Uniform, Type);                        \
        addFunction1(ZeroConstant, Type);                                      \
        addFunction1(OneConstant, Type);                                       \
        addStreamConstructableFunction1(Polynomial, Type);                     \
        addFunction1(Add, Type);                                               \
        addFunction1(Repeat, Type);                                            \
        addFunction1(Scale, Type);                                             \
        addFunction1(Shift, Type);                                             \
        addFunction1(Sine, Type);                                              \
        addFunction1(Square, Type);                                            \
        addStreamConstructableFunction1(Table, Type);                          \
        addNamedFunction1(Table, Type, tableFile);                             \
        addFunction1(UniformTable, Type);                                      \
        addFunction1(NonUniformTable, Type);                                   \
        addFunction1(Coded, Type);                                             \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
