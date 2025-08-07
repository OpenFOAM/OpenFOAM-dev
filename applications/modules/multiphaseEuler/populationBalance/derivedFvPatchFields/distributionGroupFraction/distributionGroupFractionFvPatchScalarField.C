/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "distributionGroupFractionFvPatchScalarField.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionGroupFractionFvPatchScalarField::
distributionGroupFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    groupPropertyFvScalarField(iF),
    distribution_
    (
        distribution::New(dimLength, dict.subDict("distribution"), 3, -1)
    ),
    etaPtr_(nullptr)
{}


Foam::distributionGroupFractionFvPatchScalarField::
distributionGroupFractionFvPatchScalarField
(
    const distributionGroupFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    groupPropertyFvScalarField(iF),
    distribution_(ptf.distribution_, false),
    etaPtr_(nullptr)
{}


Foam::distributionGroupFractionFvPatchScalarField::
distributionGroupFractionFvPatchScalarField
(
    const distributionGroupFractionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    groupPropertyFvScalarField(iF),
    distribution_(ptf.distribution_, false),
    etaPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distributionGroupFractionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!etaPtr_.valid())
    {
        etaPtr_.set
        (
            new scalar
            (
                popBal().etaV(i(), distribution_()).value()
            )
        );
    }

    fixedValueFvPatchScalarField::operator==(etaPtr_());

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::distributionGroupFractionFvPatchScalarField::write
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
        distributionGroupFractionFvPatchScalarField
    );
}

// ************************************************************************* //
