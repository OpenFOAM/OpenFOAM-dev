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
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

#include "AtmosphericBoundaryLayerVelocity_DimensionedFieldFunction.H"
#include "AtmosphericBoundaryLayerTurbulentKineticEnergy_DimensionedFieldFunction.H"
#include "AtmosphericBoundaryLayerTurbulentEpsilon_DimensionedFieldFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDimensionedFieldFunctions(meshType, MeshType)                      \
                                                                               \
    typedef DimensionedField<scalar, meshType, Field> scalar##MeshType;        \
    typedef DimensionedField<vector, meshType, Field> vector##MeshType;        \
                                                                               \
    namespace DimensionedFieldFunctions                                        \
    {                                                                          \
        addDimensionedFieldFunction                                            \
        (                                                                      \
            AtmosphericBoundaryLayerVelocity,                                  \
            vector##MeshType                                                   \
        );                                                                     \
    }                                                                          \
                                                                               \
    namespace DimensionedFieldFunctions                                        \
    {                                                                          \
        addDimensionedFieldFunction                                            \
        (                                                                      \
            AtmosphericBoundaryLayerTurbulentKineticEnergy,                    \
            scalar##MeshType                                                   \
        );                                                                     \
    }                                                                          \
                                                                               \
    namespace DimensionedFieldFunctions                                        \
    {                                                                          \
        addDimensionedFieldFunction                                            \
        (                                                                      \
            AtmosphericBoundaryLayerTurbulentEpsilon,                          \
            scalar##MeshType                                                   \
        );                                                                     \
    }

namespace Foam
{
    makeDimensionedFieldFunctions(fvMesh, FvMesh);
    makeDimensionedFieldFunctions(fvPatch, FvPatch);
}


// ************************************************************************* //
