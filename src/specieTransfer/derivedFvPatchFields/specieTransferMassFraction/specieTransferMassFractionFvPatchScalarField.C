/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2022 OpenFOAM Foundation
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

#include "specieTransferMassFractionFvPatchScalarField.H"
#include "specieTransferVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermophysicalTransportModel.H"
#include "fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum
    <
        specieTransferMassFractionFvPatchScalarField::property,
        4
    >::names[] =
    {
        "massFraction",
        "moleFraction",
        "molarConcentration",
        "partialPressure"
    };
}

const Foam::NamedEnum
<
    Foam::specieTransferMassFractionFvPatchScalarField::property,
    4
> Foam::specieTransferMassFractionFvPatchScalarField::propertyNames_;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::basicSpecieMixture&
Foam::specieTransferMassFractionFvPatchScalarField::composition
(
    const objectRegistry& db
)
{
    const word& name = physicalProperties::typeName;

    return db.lookupObject<fluidMulticomponentThermo>(name).composition();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specieTransferMassFractionFvPatchScalarField::
specieTransferMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_("phi"),
    UName_("U"),
    phiYp_(p.size(), 0),
    timeIndex_(-1),
    c_(0),
    property_(massFraction)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::specieTransferMassFractionFvPatchScalarField::
specieTransferMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiYp_(p.size(), 0),
    timeIndex_(-1),
    c_(dict.lookupOrDefault<scalar>("c", scalar(0))),
    property_
    (
        c_ == scalar(0)
      ? massFraction
      : propertyNames_.read(dict.lookup("property"))
    )
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::specieTransferMassFractionFvPatchScalarField::
specieTransferMassFractionFvPatchScalarField
(
    const specieTransferMassFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_),
    phiYp_(p.size(), 0),
    timeIndex_(-1),
    c_(ptf.c_),
    property_(ptf.property_)
{}


Foam::specieTransferMassFractionFvPatchScalarField::
specieTransferMassFractionFvPatchScalarField
(
    const specieTransferMassFractionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_),
    phiYp_(ptf.size(), 0),
    timeIndex_(-1),
    c_(ptf.c_),
    property_(ptf.property_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::specieTransferMassFractionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    m(phiYp_, phiYp_);
}


void Foam::specieTransferMassFractionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const specieTransferMassFractionFvPatchScalarField& tiptf =
        refCast<const specieTransferMassFractionFvPatchScalarField>(ptf);

    phiYp_.rmap(tiptf.phiYp_, addr);
}


void Foam::specieTransferMassFractionFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);

    const specieTransferMassFractionFvPatchScalarField& tiptf =
        refCast<const specieTransferMassFractionFvPatchScalarField>(ptf);

    phiYp_.reset(tiptf.phiYp_);
}


const Foam::scalarField&
Foam::specieTransferMassFractionFvPatchScalarField::phiYp() const
{
    if (timeIndex_ != this->db().time().timeIndex())
    {
        timeIndex_ = this->db().time().timeIndex();

        phiYp_ = calcPhiYp();
    }

    return phiYp_;
}


void Foam::specieTransferMassFractionFvPatchScalarField::updateCoeffs()
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

    // Set the gradient and value so that the transport and diffusion combined
    // result in the desired specie flux
    valueFraction() = phip/(phip - patch().deltaCoeffs()*AAlphaEffp);
    refValue() = *this;
    refGrad() = phip*(*this - phiYp()/uPhip)/AAlphaEffp;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::specieTransferMassFractionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<scalar>(os, "c", scalar(0), c_);
    writeEntry(os, "property", propertyNames_[property_]);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
}


// ************************************************************************* //
