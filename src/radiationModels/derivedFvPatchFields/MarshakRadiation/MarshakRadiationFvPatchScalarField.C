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

#include "MarshakRadiationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "radiationModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MarshakRadiationFvPatchScalarField::MarshakRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    radiationCoupledBase(p, dict),
    TName_(dict.lookupOrDefault<word>("T", "T"))
{
    if (dict.found("value"))
    {
        refValue() = scalarField("value", dict, p.size());
    }
    else
    {
        refValue() = 0.0;
    }

    // zero gradient
    refGrad() = 0.0;

    valueFraction() = 1.0;

    fvPatchScalarField::operator=(refValue());
}


Foam::MarshakRadiationFvPatchScalarField::MarshakRadiationFvPatchScalarField
(
    const MarshakRadiationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    radiationCoupledBase
    (
        p,
        ptf.emissivityMethod(),
        ptf.emissivity_,
        mapper
    ),
    TName_(ptf.TName_)
{}


Foam::MarshakRadiationFvPatchScalarField::MarshakRadiationFvPatchScalarField
(
    const MarshakRadiationFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    radiationCoupledBase
    (
        ptf.patch(),
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    TName_(ptf.TName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MarshakRadiationFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::map(ptf, mapper);
    radiationCoupledBase::map(ptf, mapper);
}


void Foam::MarshakRadiationFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);
    radiationCoupledBase::reset(ptf);
}


void Foam::MarshakRadiationFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Temperature field
    const scalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    // Re-calc reference value
    refValue() = 4.0*constant::physicoChemical::sigma.value()*pow4(Tp);

    // Diffusion coefficient - created by radiation model's ::updateCoeffs()
    const scalarField& gamma =
        patch().lookupPatchField<volScalarField, scalar>("gammaRad");

    const scalarField emissivity(this->emissivity());

    const scalarField Ep(emissivity/(2.0*(2.0 - emissivity)));

    // Set value fraction
    valueFraction() = 1.0/(1.0 + gamma*patch().deltaCoeffs()/Ep);

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::MarshakRadiationFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    radiationCoupledBase::write(os);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        MarshakRadiationFvPatchScalarField
    );
}


// ************************************************************************* //
