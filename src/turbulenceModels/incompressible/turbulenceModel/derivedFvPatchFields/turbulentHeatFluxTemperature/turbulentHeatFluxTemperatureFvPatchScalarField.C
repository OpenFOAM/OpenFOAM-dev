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

#include "turbulentHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // declare specialization within 'Foam' namespace
    template<>
    const char* NamedEnum
    <
        Foam::incompressible::
        turbulentHeatFluxTemperatureFvPatchScalarField::heatSourceType,
        2
    >::names[] =
    {
        "power",
        "flux"
    };
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const NamedEnum
<
    turbulentHeatFluxTemperatureFvPatchScalarField::heatSourceType,
    2
> turbulentHeatFluxTemperatureFvPatchScalarField::heatSourceTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    heatSource_(hsPower),
    q_(p.size(), 0.0),
    alphaEffName_("undefinedAlphaEff")
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    heatSource_(ptf.heatSource_),
    q_(ptf.q_, mapper),
    alphaEffName_(ptf.alphaEffName_)
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    heatSource_(heatSourceTypeNames_.read(dict.lookup("heatSource"))),
    q_("q", dict, p.size()),
    alphaEffName_(dict.lookup("alphaEff"))
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=(Field<scalar>("value", dict, p.size()));
        gradient() = Field<scalar>("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& thftpsf
)
:
    fixedGradientFvPatchScalarField(thftpsf),
    heatSource_(thftpsf.heatSource_),
    q_(thftpsf.q_),
    alphaEffName_(thftpsf.alphaEffName_)
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(thftpsf, iF),
    heatSource_(thftpsf.heatSource_),
    q_(thftpsf.q_),
    alphaEffName_(thftpsf.alphaEffName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    q_.autoMap(m);
}


void turbulentHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const turbulentHeatFluxTemperatureFvPatchScalarField& thftptf =
        refCast<const turbulentHeatFluxTemperatureFvPatchScalarField>
        (
            ptf
        );

    q_.rmap(thftptf.q_, addr);
}


void turbulentHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& alphaEffp =
        patch().lookupPatchField<volScalarField, scalar>(alphaEffName_);

    // retrieve (constant) specific heat capacity from transport dictionary
    const IOdictionary& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");
    const scalar Cp0(readScalar(transportProperties.lookup("Cp0")));

    switch (heatSource_)
    {
        case hsPower:
        {
            const scalar Ap = gSum(patch().magSf());
            gradient() = q_/(Ap*Cp0*alphaEffp);
            break;
        }
        case hsFlux:
        {
            gradient() = q_/(Cp0*alphaEffp);
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "turbulentHeatFluxTemperatureFvPatchScalarField"
                "("
                    "const fvPatch&, "
                    "const DimensionedField<scalar, volMesh>&, "
                    "const dictionary&"
                ")"
            )   << "Unknown heat source type. Valid types are: "
                << heatSourceTypeNames_ << nl << exit(FatalError);
        }
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void turbulentHeatFluxTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("heatSource") << heatSourceTypeNames_[heatSource_]
        << token::END_STATEMENT << nl;
    q_.writeEntry("q", os);
    os.writeKeyword("alphaEff") << alphaEffName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentHeatFluxTemperatureFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
