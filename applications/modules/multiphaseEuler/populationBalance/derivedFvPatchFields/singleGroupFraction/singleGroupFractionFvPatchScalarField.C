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

#include "singleGroupFractionFvPatchScalarField.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleGroupFractionFvPatchScalarField::
singleGroupFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    groupPropertyFvScalarField(iF),
    index_(dict.lookup<label>("index"))
{}


Foam::singleGroupFractionFvPatchScalarField::
singleGroupFractionFvPatchScalarField
(
    const singleGroupFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    groupPropertyFvScalarField(iF),
    index_(ptf.index_)
{}


Foam::singleGroupFractionFvPatchScalarField::
singleGroupFractionFvPatchScalarField
(
    const singleGroupFractionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    groupPropertyFvScalarField(iF),
    index_(ptf.index_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::singleGroupFractionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const populationBalanceModel& popBal = this->popBal();
    const label i = this->i();
    const label phaseiFirst = popBal.diameters()[i].iFirst();
    const label phaseiLast = popBal.diameters()[i].iLast();

    // Check the index
    if (index_ < phaseiFirst || index_ > phaseiLast)
    {
        FatalErrorInFunction
            << "Group index " << index_ << " is out of range of the group "
            << "indices associated with phase " << internalField().group()
            << " (" << phaseiFirst << " -> " << phaseiLast << ")"
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::operator==(i == index_);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::singleGroupFractionFvPatchScalarField::write
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
        singleGroupFractionFvPatchScalarField
    );
}

// ************************************************************************* //
