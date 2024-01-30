/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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

#include "semiPermeableBaffleMassFractionFvPatchScalarField.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermophysicalTransportModel.H"
#include "fluidMulticomponentThermo.H"
#include "mappedFvPatchBaseBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    specieTransferMassFractionFvPatchScalarField(p, iF, dict)
{
    mappedPatchBaseBase::validateMapForField
    (
        *this,
        iF,
        dict,
        mappedPatchBaseBase::from::differentPatch
    );
}


Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const semiPermeableBaffleMassFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    specieTransferMassFractionFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const semiPermeableBaffleMassFractionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    specieTransferMassFractionFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::semiPermeableBaffleMassFractionFvPatchScalarField::calcPhiYp() const
{
    if (c_ == scalar(0))
    {
        return tmp<scalarField>(new scalarField(patch().size(), Zero));
    }

    const word& YName = internalField().name();

    // Get the mapper and the neighbouring mesh and patch
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const fvPatch& nbrPatch = mapper.nbrFvPatch();

    const fluidThermo& thermo =
        db().lookupObject<fluidThermo>(physicalProperties::typeName);

    // Get the cell-mass fractions
    const scalarField Yc(patchInternalField());
    const scalarField nbrYc
    (
        mapper.fromNeighbour
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(YName)
           .patchInternalField()
        )
    );

    // Get the patch delta coefficients multiplied by the diffusivity
    const fluidThermophysicalTransportModel& ttm =
        db().lookupType<fluidThermophysicalTransportModel>();

    const scalarField alphaEffDeltap
    (
        ttm.kappaEff(patch().index())*patch().deltaCoeffs()
       /ttm.thermo().Cp().boundaryField()[patch().index()]
    );

    const scalarField nbrAlphaEffDeltap
    (
        mapper.fromNeighbour
        (
            ttm.kappaEff(nbrPatch.index())*nbrPatch.deltaCoeffs()
           /ttm.thermo().Cp().boundaryField()[nbrPatch.index()]
        )
    );

    // Get the specie molecular weight, if needed
    scalar Wi = NaN;
    if (property_ != massFraction)
    {
        const fluidMulticomponentThermo& thermo = this->thermo(db());
        Wi = thermo.Wi(thermo.species()[YName]).value();
    }

    // Get the mixture molecular weights, if needed
    tmp<scalarField> tW, tNbrW;
    if (property_ == moleFraction || property_ == partialPressure)
    {
        tW = thermo.W(patch().index());
        tNbrW = mapper.fromNeighbour(thermo.W(nbrPatch.index()));
    }

    // Construct coefficients that convert mass fraction to the property that
    // drives the transfer
    scalarField k(patch().size(), 1), nbrK(patch().size(), 1);
    switch(property_)
    {
        case massFraction:
            break;

        case moleFraction:
            k *= tW/Wi;
            nbrK *= tNbrW/Wi;
            break;

        case molarConcentration:
            {
                k *= thermo.rho(patch().index())/Wi;
                tmp<scalarField> nbrRhop =
                    mapper.fromNeighbour(thermo.rho(nbrPatch.index()));
                nbrK *= nbrRhop/Wi;
            }
            break;

        case partialPressure:
            {
                k *= thermo.p().boundaryField()[patch().index()]*tW/Wi;
                tmp<scalarField> nbrPp =
                    mapper.fromNeighbour
                    (
                        thermo.p().boundaryField()[nbrPatch.index()]
                    );
                nbrK *= nbrPp*tNbrW/Wi;
            }
            break;
    }

    // The transport is limited by both the difference in the amount of the
    // specie across the baffle and the rates of diffusion on either side. The
    // coefficients associated with these three transfer process are
    // harmonically averaged to represent this limiting.
    return
        patch().magSf()
       /(1/c_ + k/alphaEffDeltap + nbrK/nbrAlphaEffDeltap)
       *(k*Yc - nbrK*nbrYc);
}


void Foam::semiPermeableBaffleMassFractionFvPatchScalarField::write
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
        semiPermeableBaffleMassFractionFvPatchScalarField
    );
}

// ************************************************************************* //
