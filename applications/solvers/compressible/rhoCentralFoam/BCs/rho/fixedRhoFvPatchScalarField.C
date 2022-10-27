/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "fixedRhoFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedRhoFvPatchScalarField::fixedRhoFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    pName_("p"),
    psiName_("psi")
{}


Foam::fixedRhoFvPatchScalarField::fixedRhoFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    psiName_(dict.lookupOrDefault<word>("psi", "psi"))
{}


Foam::fixedRhoFvPatchScalarField::fixedRhoFvPatchScalarField
(
    const fixedRhoFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    psiName_(ptf.psiName_)
{}


Foam::fixedRhoFvPatchScalarField::fixedRhoFvPatchScalarField
(
    const fixedRhoFvPatchScalarField& frpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(frpsf, iF),
    pName_(frpsf.pName_),
    psiName_(frpsf.psiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedRhoFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& psip =
        patch().lookupPatchField<volScalarField, scalar>(psiName_);

    const fvPatchField<scalar>& pp =
        patch().lookupPatchField<volScalarField, scalar>(pName_);

    operator==(psip*pp);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fixedRhoFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntryIfDifferent<word>(os, "p", "p", this->pName_);
    writeEntryIfDifferent<word>(os, "psi", "psi", psiName_);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedRhoFvPatchScalarField
    );
}

// ************************************************************************* //
