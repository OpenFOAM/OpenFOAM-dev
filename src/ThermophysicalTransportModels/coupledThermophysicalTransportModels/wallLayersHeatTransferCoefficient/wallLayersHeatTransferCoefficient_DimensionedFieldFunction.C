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

#include "wallLayersHeatTransferCoefficient_DimensionedFieldFunction.H"
#include "DimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace DimensionedFieldFunctions
{
    defineTypeNameAndDebug
    (
        wallLayersHeatTransferCoefficient,
        0
    );

    typedef DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>
        scalarFvPatchDimensionedFieldFunction;
    addToRunTimeSelectionTable
    (
        scalarFvPatchDimensionedFieldFunction,
        wallLayersHeatTransferCoefficient,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DimensionedFieldFunctions::
wallLayersHeatTransferCoefficient::wallLayersHeatTransferCoefficient
(
    const dictionary& dict,
    DimensionedField<scalar, fvPatch>& field
)
:
    DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>(dict, field),
    thicknessLayers_
    (
        dict.lookup<scalarList>("thicknessLayers", dimensions::length)
    ),
    kappaLayers_
    (
        dict.lookup<scalarList>("kappaLayers", dimensions::thermalConductivity)
    )
{}


Foam::DimensionedFieldFunctions::
wallLayersHeatTransferCoefficient::wallLayersHeatTransferCoefficient
(
    const wallLayersHeatTransferCoefficient& dff,
    DimensionedField<scalar, fvPatch>& field
)
:
    DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>(dff, field),
    thicknessLayers_(dff.thicknessLayers_),
    kappaLayers_(dff.kappaLayers_)
{}


Foam::autoPtr
<
    Foam::DimensionedFieldFunction
    <
        Foam::DimensionedField<Foam::scalar, Foam::fvPatch>
    >
>
Foam::DimensionedFieldFunctions::wallLayersHeatTransferCoefficient::clone
(
    DimensionedField<scalar, fvPatch>& field
) const
{
    return autoPtr<DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>>
    (
        new wallLayersHeatTransferCoefficient
        (
            *this,
            field
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DimensionedFieldFunctions::
wallLayersHeatTransferCoefficient::evaluate()
{
    // Calculate effective heat transfer coefficient by harmonic averaging
    scalar h = 0;
    forAll(thicknessLayers_, i)
    {
        h += thicknessLayers_[i]/kappaLayers_[i];
    }
    h = 1/h;

    this->field_.primitiveFieldRef() = h;
}


void Foam::DimensionedFieldFunctions::
wallLayersHeatTransferCoefficient::write
(
    Ostream& os
) const
{
    writeEntry(os, "thicknessLayers", thicknessLayers_);
    writeEntry(os, "kappaLayers", kappaLayers_);
}


// ************************************************************************* //
