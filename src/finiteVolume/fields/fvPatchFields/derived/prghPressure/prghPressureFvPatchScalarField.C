/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "prghPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prghPressureFvPatchScalarField::
prghPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    rhoName_("rho"),
    p_(p.size(), 0.0)
{}


Foam::prghPressureFvPatchScalarField::
prghPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    rhoName_(dict.lookupOrDefault<word>("rhoName", "rho")),
    p_("p", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(p_);
    }
}


Foam::prghPressureFvPatchScalarField::
prghPressureFvPatchScalarField
(
    const prghPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_),
    p_(ptf.p_, mapper)
{}


Foam::prghPressureFvPatchScalarField::
prghPressureFvPatchScalarField
(
    const prghPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    rhoName_(ptf.rhoName_),
    p_(ptf.p_)
{}


Foam::prghPressureFvPatchScalarField::
prghPressureFvPatchScalarField
(
    const prghPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    rhoName_(ptf.rhoName_),
    p_(ptf.p_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::prghPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    p_.autoMap(m);
}


void Foam::prghPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const prghPressureFvPatchScalarField& tiptf =
        refCast<const prghPressureFvPatchScalarField>(ptf);

    p_.rmap(tiptf.p_, addr);
}


void Foam::prghPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& rhop = patch().lookupPatchField<volScalarField, scalar>
    (
        rhoName_
    );

    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    operator==(p_ - rhop*((g.value() & patch().Cf())));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::prghPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    if (rhoName_ != "rho")
    {
        os.writeKeyword("rhoName")
            << rhoName_ << token::END_STATEMENT << nl;
    }
    p_.writeEntry("p", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        prghPressureFvPatchScalarField
    );
}

// ************************************************************************* //
