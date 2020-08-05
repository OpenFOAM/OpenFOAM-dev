/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "adsorptionMassFractionFvPatchScalarField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "thermophysicalTransportModel.H"
#include "basicSpecieMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adsorptionMassFractionFvPatchScalarField::
adsorptionMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    specieTransferMassFractionFvPatchScalarField(p, iF)
{}


Foam::adsorptionMassFractionFvPatchScalarField::
adsorptionMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    specieTransferMassFractionFvPatchScalarField(p, iF, dict)
{}


Foam::adsorptionMassFractionFvPatchScalarField::
adsorptionMassFractionFvPatchScalarField
(
    const adsorptionMassFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    specieTransferMassFractionFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::adsorptionMassFractionFvPatchScalarField::
adsorptionMassFractionFvPatchScalarField
(
    const adsorptionMassFractionFvPatchScalarField& ptf
)
:
    specieTransferMassFractionFvPatchScalarField(ptf)
{}


Foam::adsorptionMassFractionFvPatchScalarField::
adsorptionMassFractionFvPatchScalarField
(
    const adsorptionMassFractionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    specieTransferMassFractionFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::adsorptionMassFractionFvPatchScalarField::calcPhiYp() const
{
    if (c_ == scalar(0))
    {
        return tmp<scalarField>(new scalarField(patch().size(), Zero));
    }

    const word& YName = internalField().name();

    const fluidThermo& thermo =
        db().lookupObject<fluidThermo>(fluidThermo::dictName);

    // Get the cell-mass fraction
    const scalarField Yc(patchInternalField());

    // Get the patch delta coefficients multiplied by the diffusivity
    const thermophysicalTransportModel& ttm =
        db().lookupObject<thermophysicalTransportModel>
        (
            thermophysicalTransportModel::typeName
        );
    const scalarField alphaEffDeltap
    (
        ttm.alphaEff(patch().index())*patch().deltaCoeffs()
    );

    // Get the specie molecular weight, if needed
    scalar Wi = NaN;
    if (property_ != massFraction)
    {
        const basicSpecieMixture& mixture = composition(db());
        Wi = mixture.Wi(mixture.species()[YName]);
    }

    // Get the mixture molecular weights, if needed
    tmp<scalarField> tW;
    if (property_ == moleFraction || property_ == partialPressure)
    {
        tW = thermo.W(patch().index());
    }

    // Construct coefficients that convert mass fraction to the property that
    // drives the transfer
    scalarField k(patch().size(), 1);
    switch(property_)
    {
        case massFraction:
            break;

        case moleFraction:
            k *= tW/Wi;
            break;

        case molarConcentration:
            k *= thermo.rho(patch().index())/Wi;
            break;

        case partialPressure:
            k *= thermo.p().boundaryField()[patch().index()]*tW/Wi;
            break;
    }

    // The transport is limited by both the difference in the amount of the
    // specie available and the rate of diffusion to the wall. The coefficients
    // associated with these two transfer process are harmonically averaged to
    // represent this limiting.
    return
        patch().magSf()
       /(1/c_ + k/alphaEffDeltap)
       *k*Yc;
}


void Foam::adsorptionMassFractionFvPatchScalarField::write
(
    Ostream& os
) const
{
    specieTransferMassFractionFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adsorptionMassFractionFvPatchScalarField
    );
}

// ************************************************************************* //
