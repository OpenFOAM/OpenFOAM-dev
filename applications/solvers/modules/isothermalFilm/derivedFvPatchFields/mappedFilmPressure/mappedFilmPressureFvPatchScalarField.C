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

#include "mappedFilmPressureFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedFilmPressureFvPatchScalarField::mappedFilmPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchField<scalar>(p, iF),
    mappedFvPatchField<scalar>(p, iF)
{}


Foam::mappedFilmPressureFvPatchScalarField::mappedFilmPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchField<scalar>(p, iF, dict),
    mappedFvPatchField<scalar>(p, iF, dict)
{}


Foam::mappedFilmPressureFvPatchScalarField::mappedFilmPressureFvPatchScalarField
(
    const mappedFilmPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchField<scalar>(ptf, p, iF, mapper),
    mappedFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::mappedFilmPressureFvPatchScalarField::mappedFilmPressureFvPatchScalarField
(
    const mappedFilmPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchField<scalar>(ptf, iF),
    mappedFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedFilmPressureFvPatchScalarField::map
(
    const fvPatchField<scalar>& ptf,
    const fvPatchFieldMapper& mapper
)
{
    zeroGradientFvPatchField<scalar>::map(ptf, mapper);
    mappedFvPatchField<scalar>::clearOut();
}


void Foam::mappedFilmPressureFvPatchScalarField::reset
(
    const fvPatchField<scalar>& ptf
)
{
    zeroGradientFvPatchField<scalar>::reset(ptf);
    mappedFvPatchField<scalar>::clearOut();
}


void Foam::mappedFilmPressureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Map the neighbouring fluid patch pressure field to this patch
    this->operator==(this->mappedValues(this->nbrPatchField()));

    // Map the patch pressure to the internal field
    UIndirectList<scalar>
    (
        const_cast<Field<scalar>&>(this->primitiveField()),
        this->patch().faceCells()
    ) = *this;

    zeroGradientFvPatchField<scalar>::updateCoeffs();
}


void Foam::mappedFilmPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    mappedFvPatchField<scalar>::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mappedFilmPressureFvPatchScalarField
    );
}


// ************************************************************************* //
