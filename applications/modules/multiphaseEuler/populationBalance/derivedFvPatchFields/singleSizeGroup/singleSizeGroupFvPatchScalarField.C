/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "singleSizeGroupFvPatchScalarField.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleSizeGroupFvPatchScalarField::singleSizeGroupFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    index_(dict.lookup<label>("index"))
{}


Foam::singleSizeGroupFvPatchScalarField::singleSizeGroupFvPatchScalarField
(
    const singleSizeGroupFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    index_(ptf.index_)
{}


Foam::singleSizeGroupFvPatchScalarField::singleSizeGroupFvPatchScalarField
(
    const singleSizeGroupFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    index_(ptf.index_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::singleSizeGroupFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const diameterModels::sizeGroup& fi =
        refCast<const diameterModels::sizeGroup>(internalField());

    // Check the index
    const label firstIndex = fi.group().sizeGroups().first().i();
    const label lastIndex = fi.group().sizeGroups().last().i();
    if (index_ < firstIndex || index_ > lastIndex)
    {
        FatalErrorInFunction
            << "Size-group index " << index_ << " is out of range of the "
            << "indices associated with phase " << internalField().group()
            << " (" << firstIndex << " -> " << lastIndex << ")"
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::operator==(fi.i() == index_);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::singleSizeGroupFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "index", index_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        singleSizeGroupFvPatchScalarField
    );
}

// ************************************************************************* //
