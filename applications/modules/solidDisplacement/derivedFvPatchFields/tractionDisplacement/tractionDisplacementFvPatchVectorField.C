/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "tractionDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), Zero),
    pressure_()
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = Zero;
}


Foam::tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_("traction", dimPressure, dict, p.size()),
    pressure_
    (
        Function1<scalar>::New
        (
            "pressure",
            db().time().userUnits(),
            dimPressure,
            dict
        )
    )
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = Zero;
}


Foam::tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(mapper(tdpvf.traction_)),
    pressure_(tdpvf.pressure_, false)
{}


Foam::tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tractionDisplacementFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
    fixedGradientFvPatchVectorField::map(ptf, mapper);

    const tractionDisplacementFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementFvPatchVectorField>(ptf);

    mapper(traction_, dmptf.traction_);
}


void Foam::tractionDisplacementFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    fixedGradientFvPatchVectorField::reset(ptf);

    const tractionDisplacementFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementFvPatchVectorField>(ptf);

    traction_.reset(dmptf.traction_);
}


void Foam::tractionDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    this->updateCoeffs(pressure_->value(db().time().value()));
}


void Foam::tractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "traction", traction_);
    writeEntry(os, db().time().userUnits(), dimPressure, pressure_());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        tractionDisplacementFvPatchVectorField
    );
}


// ************************************************************************* //
