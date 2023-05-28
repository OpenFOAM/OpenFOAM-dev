/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "filmContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::filmContactAngleFvPatchScalarField::filmContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict),
    contactAngle_(contactAngleModel::New(dict))
{}


Foam::filmContactAngleFvPatchScalarField::filmContactAngleFvPatchScalarField
(
    const filmContactAngleFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(psf, p, iF, mapper),
    contactAngle_(psf.contactAngle_, false)
{}


Foam::filmContactAngleFvPatchScalarField::filmContactAngleFvPatchScalarField
(
    const filmContactAngleFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(psf, iF),
    contactAngle_(psf.contactAngle_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::filmContactAngleFvPatchScalarField::cosTheta
(
    const fvPatchVectorField& Up,
    const vectorField& nHat
) const
{
    return contactAngle_->cosTheta(Up, nHat);
}


void Foam::filmContactAngleFvPatchScalarField::write
(
    Ostream& os
) const
{
    zeroGradientFvPatchScalarField::write(os);
    writeEntry(os, contactAngle_());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        filmContactAngleFvPatchScalarField
    );
}


// ************************************************************************* //
