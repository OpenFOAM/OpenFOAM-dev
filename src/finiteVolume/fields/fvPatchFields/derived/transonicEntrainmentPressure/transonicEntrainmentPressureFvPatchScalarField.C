/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "transonicEntrainmentPressureFvPatchScalarField.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::transonicEntrainmentPressureFvPatchScalarField::
transonicEntrainmentPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    rhoName_("rho"),
    psiName_("psi"),
    phiName_("phi"),
    gamma_(0),
    Mb_(0),
    p0_(p.size(), 0)
{}


Foam::transonicEntrainmentPressureFvPatchScalarField::
transonicEntrainmentPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    psiName_(dict.lookupOrDefault<word>("psi", "psi")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    gamma_(dict.lookup<scalar>("gamma")),
    Mb_(dict.lookupOrDefault<scalar>("Mb", 0.5)),
    p0_("p0", dict, p.size())
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
        fvPatchScalarField::operator=(p0_);
    }

    refValue() = p0_;
    refGrad() = Zero;
    valueFraction() = 1;
}


Foam::transonicEntrainmentPressureFvPatchScalarField::
transonicEntrainmentPressureFvPatchScalarField
(
    const transonicEntrainmentPressureFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    rhoName_(psf.rhoName_),
    psiName_(psf.psiName_),
    phiName_(psf.phiName_),
    gamma_(psf.gamma_),
    Mb_(psf.Mb_),
    p0_(mapper(psf.p0_))
{}


Foam::transonicEntrainmentPressureFvPatchScalarField::
transonicEntrainmentPressureFvPatchScalarField
(
    const transonicEntrainmentPressureFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    rhoName_(psf.rhoName_),
    psiName_(psf.psiName_),
    phiName_(psf.phiName_),
    gamma_(psf.gamma_),
    Mb_(psf.Mb_),
    p0_(psf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::transonicEntrainmentPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    m(p0_, p0_);
}


void Foam::transonicEntrainmentPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(psf, addr);

    const transonicEntrainmentPressureFvPatchScalarField& toppsf =
        refCast<const transonicEntrainmentPressureFvPatchScalarField>(psf);

    p0_.rmap(toppsf.p0_, addr);
}


void Foam::transonicEntrainmentPressureFvPatchScalarField::reset
(
    const fvPatchScalarField& psf
)
{
    mixedFvPatchScalarField::reset(psf);

    const transonicEntrainmentPressureFvPatchScalarField& toppsf =
        refCast<const transonicEntrainmentPressureFvPatchScalarField>(psf);

    p0_.reset(toppsf.p0_);
}


void Foam::transonicEntrainmentPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    const fvPatchField<scalar>& psip =
        patch().lookupPatchField<volScalarField, scalar>(psiName_);

    scalarField Unp(phip/patch().magSf());

    if (phi.dimensions() == dimMassFlux)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        Unp /= rhop;
    }

    // Calculate the speed of sound and Mach number at the outlet patch
    const scalarField c(sqrt(gamma_/psip));
    const scalarField Ma(max(Unp/c, scalar(0)));

    const scalar gM1ByG = (gamma_ - 1)/gamma_;

    refValue() =
        p0_
       /pow
        (
            1 - (0.5*gM1ByG)*psip*negPart(Unp)*mag(Unp),
            1/gM1ByG
        );

    valueFraction() = 1 - min(max(Ma - Mb_, scalar(0))/(1 - Mb_), scalar(1));

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::transonicEntrainmentPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "psi", "psi", psiName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntry(os, "Mb", Mb_);
    writeEntry(os, "gamma", gamma_);
    writeEntry(os, "p0", p0_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        transonicEntrainmentPressureFvPatchScalarField
    );
}

// ************************************************************************* //
