/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "turbulentMixingLengthFrequencyInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentMixingLengthFrequencyInletFvPatchScalarField::
turbulentMixingLengthFrequencyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    mixingLength_(0.0),
    kName_("undefined-k")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

turbulentMixingLengthFrequencyInletFvPatchScalarField::
turbulentMixingLengthFrequencyInletFvPatchScalarField
(
    const turbulentMixingLengthFrequencyInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    mixingLength_(ptf.mixingLength_),
    kName_(ptf.kName_)
{}

turbulentMixingLengthFrequencyInletFvPatchScalarField::
turbulentMixingLengthFrequencyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    mixingLength_(readScalar(dict.lookup("mixingLength"))),
    kName_(dict.lookupOrDefault<word>("k", "k"))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

turbulentMixingLengthFrequencyInletFvPatchScalarField::
turbulentMixingLengthFrequencyInletFvPatchScalarField
(
    const turbulentMixingLengthFrequencyInletFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    mixingLength_(ptf.mixingLength_),
    kName_(ptf.kName_)
{}

turbulentMixingLengthFrequencyInletFvPatchScalarField::
turbulentMixingLengthFrequencyInletFvPatchScalarField
(
    const turbulentMixingLengthFrequencyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    mixingLength_(ptf.mixingLength_),
    kName_(ptf.kName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentMixingLengthFrequencyInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup Cmu corresponding to the turbulence model selected
    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalar Cmu =
        turbulence.coeffDict().lookupOrDefault<scalar>("Cmu", 0.09);

    const scalar Cmu25 = pow(Cmu, 0.25);

    const fvPatchScalarField& kp =
        patch().lookupPatchField<volScalarField, scalar>(kName_);

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    this->refValue() = sqrt(kp)/(Cmu25*mixingLength_);
    this->valueFraction() = 1.0 - pos(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void turbulentMixingLengthFrequencyInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("mixingLength")
        << mixingLength_ << token::END_STATEMENT << nl;
    os.writeKeyword("phi") << this->phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("k") << kName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentMixingLengthFrequencyInletFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
