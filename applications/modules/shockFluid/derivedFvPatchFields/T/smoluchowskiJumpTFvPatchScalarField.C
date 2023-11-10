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

#include "smoluchowskiJumpTFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "basicThermo.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoluchowskiJumpTFvPatchScalarField::smoluchowskiJumpTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    psiName_(dict.lookupOrDefault<word>("psi", "psi")),
    muName_(dict.lookupOrDefault<word>("mu", "mu")),
    accommodationCoeff_(dict.lookup<scalar>("accommodationCoeff")),
    Twall_("Twall", dict, p.size()),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 1.4))
{
    if
    (
        mag(accommodationCoeff_) < small
     || mag(accommodationCoeff_) > 2.0
    )
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "unphysical accommodationCoeff specified"
            << "(0 < accommodationCoeff <= 1)" << endl
            << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::smoluchowskiJumpTFvPatchScalarField::smoluchowskiJumpTFvPatchScalarField
(
    const smoluchowskiJumpTFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    muName_(ptf.muName_),
    accommodationCoeff_(ptf.accommodationCoeff_),
    Twall_(mapper(ptf.Twall_)),
    gamma_(ptf.gamma_)
{}


Foam::smoluchowskiJumpTFvPatchScalarField::smoluchowskiJumpTFvPatchScalarField
(
    const smoluchowskiJumpTFvPatchScalarField& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF),
    accommodationCoeff_(ptpsf.accommodationCoeff_),
    Twall_(ptpsf.Twall_),
    gamma_(ptpsf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smoluchowskiJumpTFvPatchScalarField::map
(
    const fvPatchField<scalar>& ptf,
    const fieldMapper& mapper
)
{
    mixedFvPatchField<scalar>::map(ptf, mapper);

    const smoluchowskiJumpTFvPatchScalarField& ptpsf =
        refCast<const smoluchowskiJumpTFvPatchScalarField>(ptf);

    mapper(Twall_, ptpsf.Twall_);
}


void Foam::smoluchowskiJumpTFvPatchScalarField::reset
(
    const fvPatchField<scalar>& ptf
)
{
    mixedFvPatchField<scalar>::reset(ptf);

    const smoluchowskiJumpTFvPatchScalarField& ptpsf =
        refCast<const smoluchowskiJumpTFvPatchScalarField>(ptf);

    Twall_.reset(ptpsf.Twall_);
}


void Foam::smoluchowskiJumpTFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& pmu =
        patch().lookupPatchField<volScalarField, scalar>(muName_);
    const fvPatchScalarField& prho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    const fvPatchField<scalar>& ppsi =
        patch().lookupPatchField<volScalarField, scalar>(psiName_);
    const fvPatchVectorField& pU =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    // Prandtl number reading consistent with rhoCentralFoam
    const dictionary& thermophysicalProperties =
        db().lookupObject<IOdictionary>(physicalProperties::typeName);

    dimensionedScalar Pr
    (
        "Pr",
        dimless,
        thermophysicalProperties.subDict("mixture").subDict("transport")
        .lookup("Pr")
    );

    Field<scalar> C2
    (
        pmu/prho
        *sqrt(ppsi*constant::mathematical::piByTwo)
        *2.0*gamma_/Pr.value()/(gamma_ + 1.0)
        *(2.0 - accommodationCoeff_)/accommodationCoeff_
    );

    Field<scalar> aCoeff(prho.snGrad() - prho/C2);
    Field<scalar> KEbyRho(0.5*magSqr(pU));

    valueFraction() = (1.0/(1.0 + patch().deltaCoeffs()*C2));
    refValue() = Twall_;
    refGrad() = 0.0;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::smoluchowskiJumpTFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "psi", "psi", psiName_);
    writeEntryIfDifferent<word>(os, "mu", "mu", muName_);

    writeEntry(os, "accommodationCoeff", accommodationCoeff_);
    writeEntry(os, "Twall", Twall_);
    writeEntry(os, "gamma", gamma_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        smoluchowskiJumpTFvPatchScalarField
    );
}


// ************************************************************************* //
