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

#include "externalWallLayersHeatTransferCoefficient_DimensionedFieldFunction.H"
#include "DimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace DimensionedFieldFunctions
{
    defineTypeNameAndDebug
    (
        externalWallLayersHeatTransferCoefficient,
        0
    );

    typedef DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>
        scalarFvPatchDimensionedFieldFunction;
    addToRunTimeSelectionTable
    (
        scalarFvPatchDimensionedFieldFunction,
        externalWallLayersHeatTransferCoefficient,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DimensionedFieldFunctions::
externalWallLayersHeatTransferCoefficient::
externalWallLayersHeatTransferCoefficient
(
    const dictionary& dict,
    DimensionedField<scalar, fvPatch>& field
)
:
    wallLayersHeatTransferCoefficient(dict, field),
    h_
    (
        dict.lookup<scalar>
        (
            "h",
            dimensions::heatFluxDensity/dimensions::temperature
        )
    )
{}


Foam::DimensionedFieldFunctions::
externalWallLayersHeatTransferCoefficient::
externalWallLayersHeatTransferCoefficient
(
    const externalWallLayersHeatTransferCoefficient& dff,
    DimensionedField<scalar, fvPatch>& field
)
:
    wallLayersHeatTransferCoefficient(dff, field),
    h_(dff.h_)
{}


Foam::autoPtr
<
    Foam::DimensionedFieldFunction
    <
        Foam::DimensionedField<Foam::scalar, Foam::fvPatch>
    >
>
Foam::DimensionedFieldFunctions::externalWallLayersHeatTransferCoefficient::
clone
(
    DimensionedField<scalar, fvPatch>& field
) const
{
    return autoPtr<DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>>
    (
        new externalWallLayersHeatTransferCoefficient
        (
            *this,
            field
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DimensionedFieldFunctions::
externalWallLayersHeatTransferCoefficient::evaluate()
{
    wallLayersHeatTransferCoefficient::evaluate();
    this->field_.primitiveFieldRef() =
        1/(1/h_ + 1/this->field_.primitiveField());
}


void Foam::DimensionedFieldFunctions::
externalWallLayersHeatTransferCoefficient::write
(
    Ostream& os
) const
{
    wallLayersHeatTransferCoefficient::write(os);
    writeEntry(os, "h", h_);
}


// ************************************************************************* //
