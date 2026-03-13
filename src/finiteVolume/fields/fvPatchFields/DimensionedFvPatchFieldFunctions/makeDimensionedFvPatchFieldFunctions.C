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
#include "DimensionedFvPatchFields.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

#include "Zonal_DimensionedFvPatchFieldFunction.H"
#include "Function1_DimensionedFieldFunction.H"
#include "Coded_DimensionedFieldFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDimensionedFvPatchFieldFunctions(Type, nullArg)                    \
                                                                               \
    typedef DimensionedField<Type, fvPatch, Field>                             \
        DimensionedField##Type##fvPatch##Field;                                \
    defineDimensionedFieldFunction                                             \
    (                                                                          \
        DimensionedField##Type##fvPatch##Field                                 \
    );                                                                         \
    namespace DimensionedFieldFunctions                                        \
    {                                                                          \
        addDimensionedFieldFunction                                            \
        (                                                                      \
            Zonal,                                                             \
            DimensionedField##Type##fvPatch##Field                             \
        );                                                                     \
        addDimensionedFieldFunction                                            \
        (                                                                      \
            Function1,                                                         \
            DimensionedField##Type##fvPatch##Field                             \
        );                                                                     \
        addDimensionedFieldFunction                                            \
        (                                                                      \
            Coded,                                                             \
            DimensionedField##Type##fvPatch##Field                             \
        );                                                                     \
    }

namespace Foam
{
    FOR_ALL_FIELD_TYPES(makeDimensionedFvPatchFieldFunctions);
}


// ************************************************************************* //
