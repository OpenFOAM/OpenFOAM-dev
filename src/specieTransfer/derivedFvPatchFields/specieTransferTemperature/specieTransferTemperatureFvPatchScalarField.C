/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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

#include "specieTransferTemperatureFvPatchScalarField.H"
#include "specieTransferMassFractionFvPatchScalarField.H"
#include "specieTransferVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermophysicalTransportModel.H"
#include "fluidMulticomponentThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specieTransferTemperatureFvPatchScalarField::
specieTransferTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool readValue
)
:
    mixedEnergyCalculatedTemperatureFvPatchScalarField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{
    if (readValue)
    {
        fvPatchScalarField::operator==(scalarField("value", dict, p.size()));
    }
}


Foam::specieTransferTemperatureFvPatchScalarField::
specieTransferTemperatureFvPatchScalarField
(
    const specieTransferTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedEnergyCalculatedTemperatureFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_)
{}


Foam::specieTransferTemperatureFvPatchScalarField::
specieTransferTemperatureFvPatchScalarField
(
    const specieTransferTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedEnergyCalculatedTemperatureFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::tmp<Foam::scalarField>
Foam::specieTransferTemperatureFvPatchScalarField::phiHep() const
{
    typedef specieTransferMassFractionFvPatchScalarField YBCType;
    const fluidMulticomponentThermo& thermo = YBCType::thermo(db());
    const PtrList<volScalarField>& Y = thermo.Y();

    const fvPatchScalarField& Tp = *this;
    const fvPatchScalarField& pp = thermo.p().boundaryField()[patch().index()];

    // Sum up the phiHep from all the species
    tmp<scalarField> tPhiHep(new scalarField(this->size(), 0));
    scalarField& phiHep = tPhiHep.ref();
    forAll(Y, i)
    {
        const fvPatchScalarField& Yp = Y[i].boundaryField()[patch().index()];

        if (!isA<YBCType>(Yp))
        {
            FatalErrorInFunction
                << "The mass-fraction condition on patch " << patch().name()
                << " is not of type " << YBCType::typeName << "."
                << exit(FatalError);
        }

        phiHep += refCast<const YBCType>(Yp).phiYp()*thermo.hei(i, pp, Tp);
    }

    return tPhiHep;
}


void Foam::specieTransferTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the fluxes
    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);
    tmp<scalarField> uPhip =
        refCast<const specieTransferVelocityFvPatchVectorField>(Up).phip();

    const fluidThermophysicalTransportModel& ttm =
        db().lookupType<fluidThermophysicalTransportModel>();

    // Get the diffusivity
    const scalarField AAlphaEffp
    (
        patch().magSf()
       *ttm.kappaEff(patch().index())
       /ttm.thermo().Cp().boundaryField()[patch().index()]
    );

    // Get the current energy to linearise around
    const fluidThermo& thermo =
        db().lookupObject<fluidThermo>(physicalProperties::typeName);
    const scalarField& hep = thermo.he().boundaryField()[patch().index()];

    heValueFraction() = phip/(phip - patch().deltaCoeffs()*AAlphaEffp);
    heRefValue() = hep;
    heRefGrad() = phip*(hep - phiHep()/uPhip)/AAlphaEffp;

    mixedEnergyCalculatedTemperatureFvPatchScalarField::updateCoeffs();
}


void Foam::specieTransferTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        specieTransferTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
