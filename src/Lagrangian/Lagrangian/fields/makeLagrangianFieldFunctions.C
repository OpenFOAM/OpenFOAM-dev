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

#include "DimensionedFieldFunction.H"
#include "LagrangianMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeLagrangianFieldFunctions(Type, nullArg)                            \
                                                                               \
    typedef DimensionedField<Type, LagrangianMesh, Field>                      \
        DimensionedField##Type##LagrangianMesh##Field;                         \
    defineDimensionedFieldFunction                                             \
    (                                                                          \
        DimensionedField##Type##LagrangianMesh##Field                          \
    );                                                                         \
                                                                               \
    typedef DimensionedField<Type, LagrangianSubMesh, Field>                   \
        DimensionedField##Type##LagrangianSubMesh##Field;                      \
    defineDimensionedFieldFunction                                             \
    (                                                                          \
        DimensionedField##Type##LagrangianSubMesh##Field                       \
    );                                                                         \
                                                                               \
    typedef DimensionedField<Type, LagrangianMesh, LagrangianField>            \
        DimensionedField##Type##LagrangianMesh##LagrangianField;               \
    defineDimensionedFieldFunction                                             \
    (                                                                          \
        DimensionedField##Type##LagrangianMesh##LagrangianField                \
    );                                                                         \
                                                                               \
    typedef DimensionedField                                                   \
    <Type, LagrangianMesh, LagrangianPrimitiveDynamicField>                    \
    DimensionedField##Type##LagrangianMesh##LagrangianPrimitiveDynamicField;   \
    defineDimensionedFieldFunction                                             \
    (                                                                          \
        DimensionedField##Type##LagrangianMesh##LagrangianPrimitiveDynamicField\
    );                                                                         \

namespace Foam
{
    FOR_ALL_FIELD_TYPES(makeLagrangianFieldFunctions);
    makeLagrangianFieldFunctions(label, );
}

// ************************************************************************* //
