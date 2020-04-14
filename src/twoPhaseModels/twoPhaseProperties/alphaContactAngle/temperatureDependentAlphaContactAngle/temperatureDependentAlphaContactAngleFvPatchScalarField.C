/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "temperatureDependentAlphaContactAngleFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::
temperatureDependentAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(p, iF),
    TName_("T"),
    theta0_()
{}


Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::
temperatureDependentAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleFvPatchScalarField(p, iF, dict),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    theta0_(Function1<scalar>::New("theta0", dict))
{
    evaluate();
}


Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::
temperatureDependentAlphaContactAngleFvPatchScalarField
(
    const temperatureDependentAlphaContactAngleFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleFvPatchScalarField(psf, p, iF, mapper),
    TName_(psf.TName_),
    theta0_(psf.theta0_, false)
{}


Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::
temperatureDependentAlphaContactAngleFvPatchScalarField
(
    const temperatureDependentAlphaContactAngleFvPatchScalarField& psf
)
:
    alphaContactAngleFvPatchScalarField(psf),
    TName_(psf.TName_),
    theta0_(psf.theta0_, false)
{}


Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::
temperatureDependentAlphaContactAngleFvPatchScalarField
(
    const temperatureDependentAlphaContactAngleFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(psf, iF),
    TName_(psf.TName_),
    theta0_(psf.theta0_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField&,
    const fvsPatchVectorField&
) const
{
    return theta0_->value
    (
        patch().lookupPatchField<volScalarField, scalar>(TName_)
    );
}


void Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::write
(
    Ostream& os
) const
{
    alphaContactAngleFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
    writeEntry(os, theta0_());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        temperatureDependentAlphaContactAngleFvPatchScalarField
    );
}

// ************************************************************************* //
