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

#include "distributionSizeGroupFvPatchScalarField.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionSizeGroupFvPatchScalarField::
distributionSizeGroupFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    distribution_
    (
        distribution::New(dimLength, dict.subDict("distribution"), 3, -1)
    ),
    etaPtr_(nullptr)
{}


Foam::distributionSizeGroupFvPatchScalarField::
distributionSizeGroupFvPatchScalarField
(
    const distributionSizeGroupFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    distribution_(ptf.distribution_, false),
    etaPtr_(nullptr)
{}


Foam::distributionSizeGroupFvPatchScalarField::
distributionSizeGroupFvPatchScalarField
(
    const distributionSizeGroupFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    distribution_(ptf.distribution_, false),
    etaPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distributionSizeGroupFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!etaPtr_.valid())
    {
        const diameterModels::sizeGroup& fi =
            refCast<const diameterModels::sizeGroup>(internalField());

        etaPtr_.set
        (
            new scalar
            (
                fi.group().popBal().etaV(fi.i(), distribution_()).value()
            )
        );
    }

    fixedValueFvPatchScalarField::operator==(etaPtr_());

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::distributionSizeGroupFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "distribution", dimLength, distribution_(), true, false);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        distributionSizeGroupFvPatchScalarField
    );
}

// ************************************************************************* //
