/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchScalarField(p, iF, dict, true),
    phiName_
    (
        cyclicPatch().owner()
      ? dict.lookupOrDefault<word>("phi", "phi")
      : word::null
    ),
    rhoName_
    (
        cyclicPatch().owner()
      ? dict.lookupOrDefault<word>("rho", "rho")
      : word::null
    ),
    D_(cyclicPatch().owner() ? dict.lookup<scalar>("D") : NaN),
    I_(cyclicPatch().owner() ? dict.lookup<scalar>("I") : NaN),
    length_(cyclicPatch().owner() ? dict.lookup<scalar>("length") : NaN),
    relaxation_
    (
        cyclicPatch().owner()
      ? dict.lookupOrDefault<scalar>("relaxation", 1)
      : NaN
    ),
    jump0_(cyclicPatch().owner() ? jump()() : scalarField(p.size()))
{}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedJumpFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_),
    relaxation_(ptf.relaxation_),
    jump0_(ptf.jump0_)
{}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_),
    relaxation_(ptf.relaxation_),
    jump0_(ptf.jump0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousBafflePressureFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (cyclicPatch().owner())
    {
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

        jumpRef() =
            -sign(Un)
            *(
                D_*turbModel.nu(patch().index())
              + I_*0.5*magUn
             )*magUn*length_;

        if (internalField().dimensions() == dimPressure)
        {
            jumpRef() *=
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);
        }

        if (relaxation_ < 1)
        {
            jumpRef() += (1 - relaxation_)*(jump0_ - jumpRef());
        }

        jump0_ = jumpRef();

        if (debug)
        {
            scalar avePressureJump = gAverage(jumpRef());
            scalar aveVelocity = gAverage(mag(Un));

            Info<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << " Average pressure drop :" << avePressureJump
                << " Average velocity :" << aveVelocity
                << endl;
        }
    }

    fixedJumpFvPatchScalarField::updateCoeffs();
}


void Foam::porousBafflePressureFvPatchField::write(Ostream& os) const
{
    fixedJumpFvPatchScalarField::write(os);

    if (cyclicPatch().owner())
    {
        writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
        writeEntry(os, "D", D_);
        writeEntry(os, "I", I_);
        writeEntry(os, "length", length_);
        writeEntryIfDifferent(os, "relaxation", scalar(1), relaxation_);
    }

    writeEntry(os, "value", *this);
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
