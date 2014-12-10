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

#include "fixedShearStressFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    tau0_(vector::zero)
{}


Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    tau0_(dict.lookupOrDefault<vector>("tau", vector::zero))
{
    fvPatchField<vector>::operator=(patchInternalField());
}


Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fixedShearStressFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    tau0_(ptf.tau0_)
{}


Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fixedShearStressFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    tau0_(ptf.tau0_)
{}


Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fixedShearStressFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    tau0_(ptf.tau0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedShearStressFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            dimensionedInternalField().group()
        )
    );

    scalarField nuEff(turbModel.nuEff()()[patch().index()]);

    const vectorField Uc(patchInternalField());

    vector tauHat = tau0_/(mag(tau0_) + ROOTVSMALL);

    const scalarField& ry = patch().deltaCoeffs();

    operator==(tauHat*(tauHat & (tau0_*(1.0/(ry*nuEff)) + Uc)));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::fixedShearStressFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    os.writeKeyword("tau") << tau0_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        fixedShearStressFvPatchVectorField
    );
}

// ************************************************************************* //
