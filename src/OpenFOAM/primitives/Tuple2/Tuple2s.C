/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "Tuple2.H"
#include "fieldTypes.H"
#include "vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineTuple2(Type1, Type2)                                             \
    template<>                                                                 \
    const char* const Tuple2<Type1, Type2>::typeName =                         \
        "Tuple2<" #Type1 "," #Type2 ">";

#define defineScalarTypeTuple2(Type, nullArg)                                  \
    defineTuple2(scalar, Type)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    FOR_ALL_FIELD_TYPES(defineScalarTypeTuple2)

    defineTuple2(vector, scalar)

    defineTuple2(scalar, vector2D)

    defineTuple2(word, scalar)
}

// ************************************************************************* //
