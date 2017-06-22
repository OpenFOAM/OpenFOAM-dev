/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "waveAlphaFvPatchScalarField.H"
#include "waveVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "levelSet.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveAlphaFvPatchScalarField::waveAlphaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    UName_("U"),
    liquid_(true),
    inletOutlet_(true)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0;
}


Foam::waveAlphaFvPatchScalarField::waveAlphaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    liquid_(dict.lookupOrDefault<bool>("liquid", true)),
    inletOutlet_(dict.lookupOrDefault<bool>("inletOutlet", true))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchScalarField::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0;
}


Foam::waveAlphaFvPatchScalarField::waveAlphaFvPatchScalarField
(
    const waveAlphaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    liquid_(ptf.liquid_),
    inletOutlet_(ptf.inletOutlet_)
{}


Foam::waveAlphaFvPatchScalarField::waveAlphaFvPatchScalarField
(
    const waveAlphaFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    liquid_(ptf.liquid_),
    inletOutlet_(ptf.inletOutlet_)
{}

Foam::waveAlphaFvPatchScalarField::waveAlphaFvPatchScalarField
(
    const waveAlphaFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    liquid_(ptf.liquid_),
    inletOutlet_(ptf.inletOutlet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveAlphaFvPatchScalarField::alpha() const
{
    const scalar t = db().time().timeOutputValue();

    const waveVelocityFvPatchVectorField& Up =
        refCast<const waveVelocityFvPatchVectorField>
        (
            patch().lookupPatchField<volVectorField, scalar>(UName_)
        );
    const waveSuperposition& waves = Up.waves();

    return
        levelSetFraction
        (
            patch(),
            waves.height(t, patch().Cf()),
            waves.height(t, patch().patch().localPoints()),
            !liquid_
        );
}


void Foam::waveAlphaFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    refValue() = alpha();

    if (inletOutlet_)
    {
        const scalarField& phip =
            patch().lookupPatchField<surfaceScalarField, scalar>("phi");

        valueFraction() = 1 - pos0(phip);
    }
    else
    {
        valueFraction() = 1;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::waveAlphaFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<bool>(os, "inletOutlet", true, inletOutlet_);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(fvPatchScalarField, waveAlphaFvPatchScalarField);
}

// ************************************************************************* //
