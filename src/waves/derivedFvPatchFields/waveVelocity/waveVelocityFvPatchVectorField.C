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
    phiName_("phi"),
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
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
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
    phiName_(ptf.phiName_),
    waves_(ptf.waves_)
{}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const waveVelocityFvPatchVectorField& ptf
)
:
    directionMixedFvPatchVectorField(ptf),
    phiName_(ptf.phiName_),
    waves_(ptf.waves_)
{}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const waveVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    phiName_(ptf.phiName_),
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

    const vectorField UWave(U());

    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const scalarField out(pos0(phip));

    // Where inflow, fix all velocity components to values specified by the
    // wave model.
    refValue() = (1 - out)*UWave;
    valueFraction() = (1 - out)*symmTensor::I;

    // Where outflow, set the normal component of the velocity to a value
    // consistent with phi, but scale it to get the volumentic flow rate
    // specified by the wave model. Tangential components are extrapolated.
    const scalar QPhip = gSum(out*phip);
    const scalar QWave = gSum(out*(UWave & patch().Sf()));
    const vectorField nBySf(patch().Sf()/sqr(patch().magSf()));
    if (QPhip > VSMALL)
    {
        refValue() += out*(QWave/QPhip)*phip*nBySf;
    }
    else
    {
        refValue() += out*QWave*nBySf;
    }
    valueFraction() += out*sqr(patch().nf());

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}


void Foam::waveVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    directionMixedFvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
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
