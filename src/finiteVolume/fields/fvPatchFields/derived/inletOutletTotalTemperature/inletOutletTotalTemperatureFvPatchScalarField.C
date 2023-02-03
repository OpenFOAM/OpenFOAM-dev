/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "inletOutletTotalTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    UName_("U"),
    psiName_("psi"),
    gamma_(0.0),
    T0_(p.size(), 0.0)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    psiName_(dict.lookupOrDefault<word>("psi", "psi")),
    gamma_(dict.lookup<scalar>("gamma")),
    T0_("T0", dict, p.size())
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    this->refValue() = Zero;
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(T0_);
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
    const inletOutletTotalTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    T0_(mapper(ptf.T0_))
{}


Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
    const inletOutletTotalTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    T0_(tppsf.T0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inletOutletTotalTemperatureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fvPatchFieldMapper& mapper
)
{
    inletOutletFvPatchScalarField::map(ptf, mapper);

    const inletOutletTotalTemperatureFvPatchScalarField& tiptf =
        refCast<const inletOutletTotalTemperatureFvPatchScalarField>(ptf);

    mapper(T0_, tiptf.T0_);
}


void Foam::inletOutletTotalTemperatureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    inletOutletFvPatchScalarField::reset(ptf);

    const inletOutletTotalTemperatureFvPatchScalarField& tiptf =
        refCast<const inletOutletTotalTemperatureFvPatchScalarField>(ptf);

    T0_.reset(tiptf.T0_);
}


void Foam::inletOutletTotalTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    const fvPatchField<scalar>& psip =
        patch().lookupPatchField<volScalarField, scalar>(psiName_);

    scalar gM1ByG = (gamma_ - 1.0)/gamma_;

    this->refValue() =
        T0_/(1.0 + 0.5*psip*gM1ByG*neg(phip)*magSqr(Up));
    this->valueFraction() = neg(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::inletOutletTotalTemperatureFvPatchScalarField::write(Ostream& os)
const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
    writeEntryIfDifferent<word>(os, "psi", "psi", psiName_);
    writeEntry(os, "gamma", gamma_);
    writeEntry(os, "T0", T0_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        inletOutletTotalTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
