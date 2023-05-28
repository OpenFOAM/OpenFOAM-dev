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

#include "totalTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalTemperatureFvPatchScalarField::totalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    psiName_(dict.lookupOrDefault<word>("psi", "psi")),
    gamma_(dict.lookup<scalar>("gamma")),
    T0_("T0", dict, p.size())
{
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
}


Foam::totalTemperatureFvPatchScalarField::totalTemperatureFvPatchScalarField
(
    const totalTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    T0_(mapper(ptf.T0_))
{}


Foam::totalTemperatureFvPatchScalarField::totalTemperatureFvPatchScalarField
(
    const totalTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    T0_(tppsf.T0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalTemperatureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchScalarField::map(ptf, mapper);

    const totalTemperatureFvPatchScalarField& tiptf =
        refCast<const totalTemperatureFvPatchScalarField>(ptf);

    mapper(T0_, tiptf.T0_);
}


void Foam::totalTemperatureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    fixedValueFvPatchScalarField::reset(ptf);

    const totalTemperatureFvPatchScalarField& tiptf =
        refCast<const totalTemperatureFvPatchScalarField>(ptf);

    T0_.reset(tiptf.T0_);
}


void Foam::totalTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<scalar>& psip =
        patch().lookupPatchField<volScalarField, scalar>(psiName_);

    scalar gM1ByG = (gamma_ - 1.0)/gamma_;

    operator==
    (
        T0_/(1.0 + 0.5*psip*gM1ByG*neg(phip)*magSqr(Up))
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::totalTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
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
        totalTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
