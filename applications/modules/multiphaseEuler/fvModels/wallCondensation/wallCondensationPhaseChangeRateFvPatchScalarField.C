/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "wallCondensationPhaseChangeRateFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * /

Foam::wallCondensationPhaseChangeRateFvPatchScalarField::
wallCondensationPhaseChangeRateFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(p, iF),
    condensing_(p.size(), scalar(0)),
    alphatLiquid_(p.size(), scalar(0)),
    alphatVapour_(p.size(), scalar(0))
{}


Foam::wallCondensationPhaseChangeRateFvPatchScalarField::
wallCondensationPhaseChangeRateFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    calculatedFvPatchScalarField(p, iF, dict),
    condensing_("condensing", dimless, dict, p.size()),
    alphatLiquid_("alphatLiquid", dimMass/dimTime/dimLength, dict, p.size()),
    alphatVapour_("alphatVapour", dimMass/dimTime/dimLength, dict, p.size())
{}


Foam::wallCondensationPhaseChangeRateFvPatchScalarField::
wallCondensationPhaseChangeRateFvPatchScalarField
(
    const wallCondensationPhaseChangeRateFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    calculatedFvPatchScalarField(psf, p, iF, mapper),
    condensing_(mapper(psf.condensing_)),
    alphatLiquid_(mapper(psf.alphatLiquid_)),
    alphatVapour_(mapper(psf.alphatVapour_))
{}


Foam::wallCondensationPhaseChangeRateFvPatchScalarField::
wallCondensationPhaseChangeRateFvPatchScalarField
(
    const wallCondensationPhaseChangeRateFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(psf, iF),
    condensing_(psf.condensing_),
    alphatLiquid_(psf.alphatLiquid_),
    alphatVapour_(psf.alphatVapour_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalarField&
Foam::wallCondensationPhaseChangeRateFvPatchScalarField::alphatLiquid() const
{
    return alphatLiquid_;
}


const Foam::scalarField&
Foam::wallCondensationPhaseChangeRateFvPatchScalarField::alphatVapour() const
{
    return alphatVapour_;
}


void Foam::wallCondensationPhaseChangeRateFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    calculatedFvPatchScalarField::map(ptf, mapper);

    const wallCondensationPhaseChangeRateFvPatchScalarField& tiptf =
        refCast<const wallCondensationPhaseChangeRateFvPatchScalarField>(ptf);

    mapper(condensing_, tiptf.condensing_);
    mapper(alphatLiquid_, tiptf.alphatLiquid_);
    mapper(alphatVapour_, tiptf.alphatVapour_);
}


void Foam::wallCondensationPhaseChangeRateFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    calculatedFvPatchScalarField::reset(ptf);

    const wallCondensationPhaseChangeRateFvPatchScalarField& tiptf =
        refCast<const wallCondensationPhaseChangeRateFvPatchScalarField>(ptf);

    condensing_.reset(tiptf.condensing_);
    alphatLiquid_.reset(tiptf.alphatLiquid_);
    alphatVapour_.reset(tiptf.alphatVapour_);
}


void Foam::wallCondensationPhaseChangeRateFvPatchScalarField::updateCoeffs()
{
    NotImplemented;
}


void Foam::wallCondensationPhaseChangeRateFvPatchScalarField::write
(
    Ostream& os
) const
{
    calculatedFvPatchScalarField::write(os);

    writeEntry(os, "condensing", condensing_);
    writeEntry(os, "alphatLiquid", alphatLiquid_);
    writeEntry(os, "alphatVapour", alphatVapour_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeNullConstructablePatchTypeField
    (
        fvPatchScalarField,
        wallCondensationPhaseChangeRateFvPatchScalarField
    );
}


// ************************************************************************* //
