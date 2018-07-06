/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "outletPhaseMeanVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outletPhaseMeanVelocityFvPatchVectorField
::outletPhaseMeanVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(p, iF),
    Umean_(0),
    ramp_(),
    alphaName_("none")
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::outletPhaseMeanVelocityFvPatchVectorField
::outletPhaseMeanVelocityFvPatchVectorField
(
    const outletPhaseMeanVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper),
    Umean_(ptf.Umean_),
    ramp_(ptf.ramp_, false),
    alphaName_(ptf.alphaName_)
{}


Foam::outletPhaseMeanVelocityFvPatchVectorField
::outletPhaseMeanVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF),
    Umean_(readScalar(dict.lookup("Umean"))),
    ramp_
    (
        dict.found("ramp")
      ? Function1<scalar>::New("ramp", dict)
      : autoPtr<Function1<scalar>>()
    ),
    alphaName_(dict.lookup("alpha"))
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }
}


Foam::outletPhaseMeanVelocityFvPatchVectorField
::outletPhaseMeanVelocityFvPatchVectorField
(
    const outletPhaseMeanVelocityFvPatchVectorField& ptf
)
:
    mixedFvPatchField<vector>(ptf),
    Umean_(ptf.Umean_),
    ramp_(ptf.ramp_, false),
    alphaName_(ptf.alphaName_)
{}


Foam::outletPhaseMeanVelocityFvPatchVectorField
::outletPhaseMeanVelocityFvPatchVectorField
(
    const outletPhaseMeanVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF),
    Umean_(ptf.Umean_),
    ramp_(ptf.ramp_, false),
    alphaName_(ptf.alphaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::outletPhaseMeanVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();

    scalarField alphap =
        patch().lookupPatchField<volScalarField, scalar>(alphaName_);

    alphap = max(alphap, scalar(0));
    alphap = min(alphap, scalar(1));

    // Get the patchInternalField (zero-gradient field)
    vectorField Uzg(patchInternalField());

    // Calculate the phase mean zero-gradient velocity
    scalar Uzgmean =
        gSum(alphap*(patch().Sf() & Uzg))
       /gSum(alphap*patch().magSf());

    // Set the refValue and valueFraction to adjust the boundary field
    // such that the phase mean is Umean_
    const scalar Umean = (ramp_.valid() ? ramp_->value(t) : 1)*Umean_;
    if (Uzgmean >= Umean)
    {
        refValue() = Zero;
        valueFraction() = 1.0 - Umean/Uzgmean;
    }
    else
    {
        refValue() = (Umean + Uzgmean)*patch().nf();
        valueFraction() = 1.0 - Uzgmean/Umean;
    }

    mixedFvPatchField<vector>::updateCoeffs();
}


void Foam::outletPhaseMeanVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);

    os.writeKeyword("Umean") << Umean_ << token::END_STATEMENT << nl;
    if (ramp_.valid())
    {
        ramp_->writeData(os);
    }
    os.writeKeyword("alpha") << alphaName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       outletPhaseMeanVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
