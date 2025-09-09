/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "wallBoilingPhaseChangeRateFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

const hashedWordList
wallBoilingPhaseChangeRateFvPatchScalarField::propertyNames_
({
    "wetFraction",
    "dDeparture",
    "fDeparture",
    "nucleationSiteDensity",
    "qQuenching",
    "qEvaporative"
});

const List<const scalarField wallBoilingPhaseChangeRateFvPatchScalarField::*>
wallBoilingPhaseChangeRateFvPatchScalarField::propertyPtrs_
({
    &wallBoilingPhaseChangeRateFvPatchScalarField::wetFraction_,
    &wallBoilingPhaseChangeRateFvPatchScalarField::dDeparture_,
    &wallBoilingPhaseChangeRateFvPatchScalarField::fDeparture_,
    &wallBoilingPhaseChangeRateFvPatchScalarField::nucleationSiteDensity_,
    &wallBoilingPhaseChangeRateFvPatchScalarField::qQuenching_,
    &wallBoilingPhaseChangeRateFvPatchScalarField::qEvaporative_
});

const dimensionSet dimInvArea(inv(dimArea));
const dimensionSet dimHeatFlux(dimEnergy*inv(dimTime*dimArea));

const List<const dimensionSet*>
wallBoilingPhaseChangeRateFvPatchScalarField::propertyDimensions_
({
    &dimless,
    &dimLength,
    &dimRate,
    &dimInvArea,
    &dimHeatFlux,
    &dimHeatFlux
});

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * /

Foam::wallBoilingPhaseChangeRateFvPatchScalarField::
wallBoilingPhaseChangeRateFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(p, iF),
    wetFraction_(p.size(), scalar(0)),
    dDeparture_(p.size(), vGreat),
    fDeparture_(p.size(), scalar(0)),
    nucleationSiteDensity_(p.size(), scalar(0)),
    qQuenching_(p.size(), scalar(0)),
    qEvaporative_(p.size(), scalar(0)),
    alphatLiquid_(p.size(), scalar(0)),
    alphatVapour_(p.size(), scalar(0))
{}


Foam::wallBoilingPhaseChangeRateFvPatchScalarField::
wallBoilingPhaseChangeRateFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    calculatedFvPatchScalarField(p, iF, dict),
    wetFraction_("wetFraction", dimless, dict, p.size()),
    dDeparture_("dDeparture", dimLength, dict, p.size()),
    fDeparture_("fDeparture", dimRate, dict, p.size()),
    nucleationSiteDensity_("nucleationSiteDensity", dimInvArea, dict, p.size()),
    qQuenching_("qQuenching", dimHeatFlux, dict, p.size()),
    qEvaporative_("qEvaporative", dimHeatFlux, dict, p.size()),
    alphatLiquid_("alphatLiquid", dimMass/dimTime/dimLength, dict, p.size()),
    alphatVapour_("alphatVapour", dimMass/dimTime/dimLength, dict, p.size())
{}


Foam::wallBoilingPhaseChangeRateFvPatchScalarField::
wallBoilingPhaseChangeRateFvPatchScalarField
(
    const wallBoilingPhaseChangeRateFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    calculatedFvPatchScalarField(psf, p, iF, mapper),
    wetFraction_(mapper(psf.wetFraction_)),
    dDeparture_(mapper(psf.dDeparture_)),
    fDeparture_(mapper(psf.fDeparture_)),
    nucleationSiteDensity_(mapper(psf.nucleationSiteDensity_)),
    qQuenching_(mapper(psf.qQuenching_)),
    qEvaporative_(mapper(psf.qEvaporative_)),
    alphatLiquid_(mapper(psf.alphatLiquid_)),
    alphatVapour_(mapper(psf.alphatVapour_))
{}


Foam::wallBoilingPhaseChangeRateFvPatchScalarField::
wallBoilingPhaseChangeRateFvPatchScalarField
(
    const wallBoilingPhaseChangeRateFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(psf, iF),
    wetFraction_(psf.wetFraction_),
    dDeparture_(psf.dDeparture_),
    fDeparture_(psf.fDeparture_),
    nucleationSiteDensity_(psf.nucleationSiteDensity_),
    qQuenching_(psf.qQuenching_),
    qEvaporative_(psf.qEvaporative_),
    alphatLiquid_(psf.alphatLiquid_),
    alphatVapour_(psf.alphatVapour_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalarField&
Foam::wallBoilingPhaseChangeRateFvPatchScalarField::property
(
    const word& name
) const
{
    return this->*propertyPtrs_[propertyNames_[name]];
}


const Foam::dimensionSet&
Foam::wallBoilingPhaseChangeRateFvPatchScalarField::propertyDimensions
(
    const word& name
)
{
    return *propertyDimensions_[propertyNames_[name]];
}


const Foam::scalarField&
Foam::wallBoilingPhaseChangeRateFvPatchScalarField::alphatLiquid() const
{
    return alphatLiquid_;
}


const Foam::scalarField&
Foam::wallBoilingPhaseChangeRateFvPatchScalarField::alphatVapour() const
{
    return alphatVapour_;
}


void Foam::wallBoilingPhaseChangeRateFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    calculatedFvPatchScalarField::map(ptf, mapper);

    const wallBoilingPhaseChangeRateFvPatchScalarField& tiptf =
        refCast<const wallBoilingPhaseChangeRateFvPatchScalarField>(ptf);

    mapper(wetFraction_, tiptf.wetFraction_);
    mapper(dDeparture_, tiptf.dDeparture_);
    mapper(fDeparture_, tiptf.fDeparture_);
    mapper(nucleationSiteDensity_, tiptf.nucleationSiteDensity_);
    mapper(qQuenching_, tiptf.qQuenching_);
    mapper(qEvaporative_, tiptf.qEvaporative_);
    mapper(alphatLiquid_, tiptf.alphatLiquid_);
    mapper(alphatVapour_, tiptf.alphatVapour_);
}


void Foam::wallBoilingPhaseChangeRateFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    calculatedFvPatchScalarField::reset(ptf);

    const wallBoilingPhaseChangeRateFvPatchScalarField& tiptf =
        refCast<const wallBoilingPhaseChangeRateFvPatchScalarField>(ptf);

    wetFraction_.reset(tiptf.wetFraction_);
    dDeparture_.reset(tiptf.dDeparture_);
    fDeparture_.reset(tiptf.fDeparture_);
    nucleationSiteDensity_.reset(tiptf.nucleationSiteDensity_);
    qQuenching_.reset(tiptf.qQuenching_);
    qEvaporative_.reset(tiptf.qEvaporative_);
    alphatLiquid_.reset(tiptf.alphatLiquid_);
    alphatVapour_.reset(tiptf.alphatVapour_);
}


void Foam::wallBoilingPhaseChangeRateFvPatchScalarField::updateCoeffs()
{
    NotImplemented;
}


void Foam::wallBoilingPhaseChangeRateFvPatchScalarField::write
(
    Ostream& os
) const
{
    calculatedFvPatchScalarField::write(os);

    writeEntry(os, "wetFraction", wetFraction_);
    writeEntry(os, "dDeparture", dDeparture_);
    writeEntry(os, "fDeparture", fDeparture_);
    writeEntry(os, "nucleationSiteDensity", nucleationSiteDensity_);
    writeEntry(os, "qQuenching", qQuenching_);
    writeEntry(os, "qEvaporative", qEvaporative_);
    writeEntry(os, "alphatLiquid", alphatLiquid_);
    writeEntry(os, "alphatVapour", alphatVapour_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeNullConstructablePatchTypeField
    (
        fvPatchScalarField,
        wallBoilingPhaseChangeRateFvPatchScalarField
    );
}


// ************************************************************************* //
