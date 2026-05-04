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

#include "displacementGapHeatTransferCoefficient_DimensionedFieldFunction.H"
#include "DimensionedFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace DimensionedFieldFunctions
{
    defineTypeNameAndDebug
    (
        displacementGapHeatTransferCoefficient,
        0
    );

    typedef DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>
        DimensionedFieldFunctionScalarFvPatch;
    addToRunTimeSelectionTable
    (
        DimensionedFieldFunctionScalarFvPatch,
        displacementGapHeatTransferCoefficient,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DimensionedFieldFunctions::
displacementGapHeatTransferCoefficient::displacementGapHeatTransferCoefficient
(
    const dictionary& dict,
    DimensionedField<scalar, fvPatch>& field
)
:
    DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>(dict, field),
    kappa_("kappa", dimThermalConductivity, dict)
{}


Foam::DimensionedFieldFunctions::
displacementGapHeatTransferCoefficient::displacementGapHeatTransferCoefficient
(
    const displacementGapHeatTransferCoefficient& dff,
    DimensionedField<scalar, fvPatch>& field
)
:
    DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>(dff, field),
    kappa_(dff.kappa_)
{}


Foam::autoPtr
<
    Foam::DimensionedFieldFunction
    <
        Foam::DimensionedField<Foam::scalar, Foam::fvPatch>
    >
>
Foam::DimensionedFieldFunctions::displacementGapHeatTransferCoefficient::clone
(
    DimensionedField<scalar, fvPatch>& field
) const
{
    return autoPtr<DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>>
    (
        new displacementGapHeatTransferCoefficient
        (
            *this,
            field
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DimensionedFieldFunctions::
displacementGapHeatTransferCoefficient::evaluate()
{}


bool Foam::DimensionedFieldFunctions::
displacementGapHeatTransferCoefficient::update()
{
    const fvPatch& p = this->field_.mesh();

    // Lookup the displacement for the patch
    const fvPatchVectorField& Dp
        = p.lookupPatchField<volVectorField, vector>("D");

    this->field_.primitiveFieldRef() = kappa_.value()/max(-p.nf() & Dp, small);

    return true;
}


void Foam::DimensionedFieldFunctions::
displacementGapHeatTransferCoefficient::write
(
    Ostream& os
) const
{
    writeEntry(os, kappa_);
}


// ************************************************************************* //
