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
#include "pointMesh.H"
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"

#include "TimeFunction_DimensionedFieldFunction.H"
#include "Zonal_DimensionedFieldFunction.H"
#include "DistanceFunction_DimensionedFieldFunction.H"
#include "Coded_DimensionedFieldFunction.H"
#include "FieldFunction_DimensionedFieldFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDimensionedPointFieldFunctions(Type, nullArg)                      \
                                                                               \
    typedef DimensionedField<Type, pointMesh, Field> Type##PointMesh;          \
    defineDimensionedFieldFunction(Type##PointMesh);                           \
    namespace DimensionedFieldFunctions                                        \
    {                                                                          \
        addDimensionedFieldFunction(TimeFunction, Type##PointMesh);            \
        addDimensionedFieldFunction(Zonal, Type##PointMesh);                   \
        addDimensionedFieldFunction(DistanceFunction, Type##PointMesh);        \
        addDimensionedFieldFunction(Coded, Type##PointMesh);                   \
        addDimensionedFieldFunction(FieldFunction, Type##PointMesh);           \
    }

namespace Foam
{
    FOR_ALL_FIELD_TYPES(makeDimensionedPointFieldFunctions);
}


// ************************************************************************* //
