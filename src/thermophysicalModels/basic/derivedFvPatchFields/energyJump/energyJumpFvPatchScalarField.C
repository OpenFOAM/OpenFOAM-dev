/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "energyJumpFvPatchScalarField.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    jumpCyclicFvPatchScalarField(p, iF),
    jump_(p.size(), Zero)
{
    evaluateNoUpdateCoeffs();
}


Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    jumpCyclicFvPatchScalarField(p, iF, dict),
    jump_("jump", iF.dimensions(), dict, p.size())
{
    evaluateNoUpdateCoeffs();
}


Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const energyJumpFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    jumpCyclicFvPatchScalarField(ptf, p, iF, mapper),
    jump_(mapper(ptf.jump_))
{}


Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const energyJumpFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    jumpCyclicFvPatchScalarField(ptf, iF),
    jump_(ptf.jump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::energyJumpFvPatchScalarField::jump() const
{
    return jump_;
}


void Foam::energyJumpFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    jumpCyclicFvPatchScalarField::map(ptf, mapper);

    const energyJumpFvPatchScalarField& tiptf =
        refCast<const energyJumpFvPatchScalarField>(ptf);

    mapper(jump_, tiptf.jump_);
}


void Foam::energyJumpFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    jumpCyclicFvPatchScalarField::reset(ptf);

    const energyJumpFvPatchScalarField& tiptf =
        refCast<const energyJumpFvPatchScalarField>(ptf);

    jump_.reset(tiptf.jump_);
}


void Foam::energyJumpFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const basicThermo& thermo = basicThermo::lookupThermo(*this);

    jumpCyclicFvPatchScalarField& Tp =
        const_cast<jumpCyclicFvPatchScalarField&>
        (
            refCast<const jumpCyclicFvPatchScalarField>
            (
                thermo.T().boundaryField()[patch().index()]
            )
        );

    // Force update of the jump
    Tp.updateCoeffs();

    // Convert the temperature jump to an energy jump
    jump_ = thermo.he(Tp.jump(), patch().faceCells());

    jumpCyclicFvPatchScalarField::updateCoeffs();
}


void Foam::energyJumpFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "jump", jump_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makeNullConstructablePatchTypeField
   (
       fvPatchScalarField,
       energyJumpFvPatchScalarField
   );
}


// ************************************************************************* //
