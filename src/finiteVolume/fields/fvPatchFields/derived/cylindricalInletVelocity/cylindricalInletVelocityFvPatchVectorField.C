/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "cylindricalInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    centre_(pTraits<vector>::zero),
    axis_(pTraits<vector>::zero),
    axialVelocity_(),
    radialVelocity_(),
    rpm_()
{}


Foam::cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const cylindricalInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    centre_(ptf.centre_),
    axis_(ptf.axis_),
    axialVelocity_(ptf.axialVelocity_().clone().ptr()),
    radialVelocity_(ptf.radialVelocity_().clone().ptr()),
    rpm_(ptf.rpm_().clone().ptr())
{}


Foam::cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    centre_(dict.lookup("centre")),
    axis_(dict.lookup("axis")),
    axialVelocity_(DataEntry<scalar>::New("axialVelocity", dict)),
    radialVelocity_(DataEntry<scalar>::New("radialVelocity", dict)),
    rpm_(DataEntry<scalar>::New("rpm", dict))
{}


Foam::cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const cylindricalInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    centre_(ptf.centre_),
    axis_(ptf.axis_),
    axialVelocity_(ptf.axialVelocity_().clone().ptr()),
    radialVelocity_(ptf.radialVelocity_().clone().ptr()),
    rpm_(ptf.rpm_().clone().ptr())
{}


Foam::cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const cylindricalInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    centre_(ptf.centre_),
    axis_(ptf.axis_),
    axialVelocity_(ptf.axialVelocity_().clone().ptr()),
    radialVelocity_(ptf.radialVelocity_().clone().ptr()),
    rpm_(ptf.rpm_().clone().ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cylindricalInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    const scalar axialVelocity = axialVelocity_->value(t);
    const scalar radialVelocity = radialVelocity_->value(t);
    const scalar rpm = rpm_->value(t);

    vector hatAxis = axis_/mag(axis_);

    const vectorField r(patch().Cf() - centre_);
    const vectorField d(r - (hatAxis & r)*hatAxis);

    tmp<vectorField> tangVel
    (
        (rpm*constant::mathematical::pi/30.0)*(hatAxis) ^ d
    );

    operator==(tangVel + hatAxis*axialVelocity + radialVelocity*d/mag(d));

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::cylindricalInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("centre") << centre_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    axialVelocity_->writeData(os);
    radialVelocity_->writeData(os);
    rpm_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       cylindricalInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
