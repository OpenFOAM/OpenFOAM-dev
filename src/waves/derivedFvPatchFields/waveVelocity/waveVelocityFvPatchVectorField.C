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

#include "waveVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "levelSet.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    waves_(db())
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    waves_(db(), dict)
{
    if (dict.found("value"))
    {
        fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const waveVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    waves_(ptf.waves_)
{}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const waveVelocityFvPatchVectorField& ptf
)
:
    directionMixedFvPatchVectorField(ptf),
    waves_(ptf.waves_)
{}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const waveVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    waves_(ptf.waves_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::waveVelocityFvPatchVectorField::U() const
{
    const scalar t = db().time().timeOutputValue();

    return
        levelSetAverage
        (
            patch(),
            waves_.height(t, patch().Cf()),
            waves_.height(t, patch().patch().localPoints()),
            waves_.UGas(t, patch().Cf())(),
            waves_.UGas(t, patch().patch().localPoints())(),
            waves_.ULiquid(t, patch().Cf())(),
            waves_.ULiquid(t, patch().patch().localPoints())()
        );
}


void Foam::waveVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    refValue() = U();

    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    const symmTensorField nnp(sqr(patch().nf()));

    // Fix all components if inflow, just normal if outflow
    valueFraction() = (1 - pos0(phip))*I + pos0(phip)*nnp;

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}


void Foam::waveVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    directionMixedFvPatchVectorField::write(os);
    waves_.write(os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        waveVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
