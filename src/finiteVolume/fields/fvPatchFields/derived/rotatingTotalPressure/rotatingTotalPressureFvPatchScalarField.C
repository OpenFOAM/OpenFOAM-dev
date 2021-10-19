/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "rotatingTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rotatingTotalPressureFvPatchScalarField::
rotatingTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(p, iF),
    omega_()
{}


Foam::rotatingTotalPressureFvPatchScalarField::
rotatingTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    totalPressureFvPatchScalarField(p, iF, dict),
    omega_(Function1<vector>::New("omega", dict))
{}


Foam::rotatingTotalPressureFvPatchScalarField::
rotatingTotalPressureFvPatchScalarField
(
    const rotatingTotalPressureFvPatchScalarField& rtppsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    totalPressureFvPatchScalarField(rtppsf, p, iF, mapper),
    omega_(rtppsf.omega_, false)
{}


Foam::rotatingTotalPressureFvPatchScalarField::
rotatingTotalPressureFvPatchScalarField
(
    const rotatingTotalPressureFvPatchScalarField& rtppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(rtppsf, iF),
    omega_(rtppsf.omega_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rotatingTotalPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().userTimeValue();
    const vector omega = omega_->value(t);
    const vector axis = omega/mag(omega);

    const vectorField Up
    (
        patch().lookupPatchField<volVectorField, vector>(UName_)
      + (omega ^ (patch().Cf() - axis*(axis & patch().Cf())))
    );

    totalPressureFvPatchScalarField::updateCoeffs(p0_, Up);
}


void Foam::rotatingTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    totalPressureFvPatchScalarField::write(os);
    writeEntry(os, omega_());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        rotatingTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
