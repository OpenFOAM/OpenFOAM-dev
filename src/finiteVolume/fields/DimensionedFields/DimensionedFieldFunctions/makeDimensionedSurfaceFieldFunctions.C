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
#include "fvMesh.H"
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"

#include "Zonal_DimensionedFieldFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDimensionedSurfaceFieldFunctions(Type, nullArg)                    \
                                                                               \
    typedef DimensionedField<Type, surfaceMesh, Field>                         \
        DimensionedField##Type##surfaceMesh##Field;                            \
    defineDimensionedFieldFunction                                             \
    (                                                                          \
        DimensionedField##Type##surfaceMesh##Field                             \
    );                                                                         \
    namespace DimensionedFieldFunctions                                        \
    {                                                                          \
        addDimensionedFieldFunction                                            \
        (                                                                      \
            Zonal,                                                             \
            DimensionedField##Type##surfaceMesh##Field                         \
        );                                                                     \
    }

namespace Foam
{
    FOR_ALL_FIELD_TYPES(makeDimensionedSurfaceFieldFunctions);
    makeDimensionedSurfaceFieldFunctions(label, );
}

// ************************************************************************* //
