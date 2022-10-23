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

#include "porousBafflePressureFvPatchField.H"
#include "surfaceFields.H"
#include "momentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    D_(0),
    I_(0),
    length_(0),
    relaxation_(1),
    jump0_(jump_)
{}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<scalar>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    D_(dict.lookup<scalar>("D")),
    I_(dict.lookup<scalar>("I")),
    length_(dict.lookup<scalar>("length")),
    relaxation_(dict.lookupOrDefault<scalar>("relaxation", 1)),
    jump0_(jump_)
{
    fvPatchField<scalar>::operator=
    (
        Field<scalar>("value", dict, p.size())
    );
}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_),
    relaxation_(ptf.relaxation_),
    jump0_(ptf.jump_)
{}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_),
    relaxation_(ptf.relaxation_),
    jump0_(ptf.jump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousBafflePressureFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    scalarField Un(phip/patch().magSf());

    if (phi.dimensions() == dimMassFlux)
    {
        Un /= patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    }

    const scalarField magUn(mag(Un));

    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    jump_ =
        -sign(Un)
        *(
            D_*turbModel.nu(patch().index())
          + I_*0.5*magUn
         )*magUn*length_;

    if (internalField().dimensions() == dimPressure)
    {
        jump_ *= patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    }

    if (relaxation_ < 1)
    {
        jump_ += (1 - relaxation_)*(jump0_ - jump_);
    }

    jump0_ = jump_;

    if (debug)
    {
        scalar avePressureJump = gAverage(jump_);
        scalar aveVelocity = gAverage(mag(Un));

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << " Average pressure drop :" << avePressureJump
            << " Average velocity :" << aveVelocity
            << endl;
    }

    fixedJumpFvPatchField<scalar>::updateCoeffs();
}


void Foam::porousBafflePressureFvPatchField::write(Ostream& os) const
{
    fixedJumpFvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntry(os, "D", D_);
    writeEntry(os, "I", I_);
    writeEntry(os, "length", length_);
    writeEntryIfDifferent(os, "relaxation", scalar(1), relaxation_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        porousBafflePressureFvPatchField
    );
}

// ************************************************************************* //
