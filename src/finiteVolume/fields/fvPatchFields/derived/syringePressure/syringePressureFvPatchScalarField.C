/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "syringePressureFvPatchScalarField.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::syringePressureFvPatchScalarField::syringePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_("phi"),
    curTimeIndex_(-1)
{}


Foam::syringePressureFvPatchScalarField::syringePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    Ap_(dict.lookup<scalar>("Ap")),
    Sp_(dict.lookup<scalar>("Sp")),
    VsI_(dict.lookup<scalar>("VsI")),
    tas_(dict.lookup<scalar>("tas")),
    tae_(dict.lookup<scalar>("tae")),
    tds_(dict.lookup<scalar>("tds")),
    tde_(dict.lookup<scalar>("tde")),
    psI_(dict.lookup<scalar>("psI")),
    psi_(dict.lookup<scalar>("psi")),
    ams_(dict.lookup<scalar>("ams")),
    ams0_(ams_),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    curTimeIndex_(-1)
{
    scalar ps = (psI_*VsI_ + ams_/psi_)/Vs(db().time().value());
    fvPatchField<scalar>::operator=(ps);
}


Foam::syringePressureFvPatchScalarField::syringePressureFvPatchScalarField
(
    const syringePressureFvPatchScalarField& sppsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(sppsf, p, iF, mapper),
    Ap_(sppsf.Ap_),
    Sp_(sppsf.Sp_),
    VsI_(sppsf.VsI_),
    tas_(sppsf.tas_),
    tae_(sppsf.tae_),
    tds_(sppsf.tds_),
    tde_(sppsf.tde_),
    psI_(sppsf.psI_),
    psi_(sppsf.psi_),
    ams_(sppsf.ams_),
    ams0_(sppsf.ams0_),
    phiName_(sppsf.phiName_),
    curTimeIndex_(-1)
{}


Foam::syringePressureFvPatchScalarField::syringePressureFvPatchScalarField
(
    const syringePressureFvPatchScalarField& sppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(sppsf, iF),
    Ap_(sppsf.Ap_),
    Sp_(sppsf.Sp_),
    VsI_(sppsf.VsI_),
    tas_(sppsf.tas_),
    tae_(sppsf.tae_),
    tds_(sppsf.tds_),
    tde_(sppsf.tde_),
    psI_(sppsf.psI_),
    psi_(sppsf.psi_),
    ams_(sppsf.ams_),
    ams0_(sppsf.ams0_),
    phiName_(sppsf.phiName_),
    curTimeIndex_(-1)
{}


Foam::syringePressureFvPatchScalarField::syringePressureFvPatchScalarField
(
    const syringePressureFvPatchScalarField& sppsf
)
:
    fixedValueFvPatchScalarField(sppsf),
    Ap_(sppsf.Ap_),
    Sp_(sppsf.Sp_),
    VsI_(sppsf.VsI_),
    tas_(sppsf.tas_),
    tae_(sppsf.tae_),
    tds_(sppsf.tds_),
    tde_(sppsf.tde_),
    psI_(sppsf.psI_),
    psi_(sppsf.psi_),
    ams_(sppsf.ams_),
    ams0_(sppsf.ams0_),
    phiName_(sppsf.phiName_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::syringePressureFvPatchScalarField::Vs(const scalar t) const
{
    if (t < tas_)
    {
        return VsI_;
    }
    else if (t < tae_)
    {
        return
            VsI_
          + 0.5*Ap_*Sp_*sqr(t - tas_)/(tae_ - tas_);
    }
    else if (t < tds_)
    {
        return
            VsI_
          + 0.5*Ap_*Sp_*(tae_ - tas_)
          + Ap_*Sp_*(t - tae_);
    }
    else if (t < tde_)
    {
        return
            VsI_
          + 0.5*Ap_*Sp_*(tae_ - tas_)
          + Ap_*Sp_*(tds_ - tae_)
          + Ap_*Sp_*(t - tds_)
          - 0.5*Ap_*Sp_*sqr(t - tds_)/(tde_ - tds_);
    }
    else
    {
        return
            VsI_
          + 0.5*Ap_*Sp_*(tae_ - tas_)
          + Ap_*Sp_*(tds_ - tae_)
          + 0.5*Ap_*Sp_*(tde_ - tds_);
    }
}


void Foam::syringePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        ams0_ = ams_;
        curTimeIndex_ = db().time().timeIndex();
    }

    scalar t = db().time().value();
    scalar deltaT = db().time().deltaTValue();

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        ams_ = ams0_ + deltaT*sum((*this*psi_)*phip);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        ams_ = ams0_ + deltaT*sum(phip);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    scalar ps = (psI_*VsI_ + ams_/psi_)/Vs(t);

    operator==(ps);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::syringePressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "Ap", Ap_);
    writeEntry(os, "Sp", Sp_);
    writeEntry(os, "VsI", VsI_);
    writeEntry(os, "tas", tas_);
    writeEntry(os, "tae", tae_);
    writeEntry(os, "tds", tds_);
    writeEntry(os, "tde", tde_);
    writeEntry(os, "psI", psI_);
    writeEntry(os, "psi", psi_);
    writeEntry(os, "ams", ams_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        syringePressureFvPatchScalarField
    );
}

// ************************************************************************* //
