/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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
#include "mappedFvPatchBaseBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedFilmPressureFvPatchScalarField::mappedFilmPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict)
{}


Foam::mappedFilmPressureFvPatchScalarField::mappedFilmPressureFvPatchScalarField
(
    const mappedFilmPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::mappedFilmPressureFvPatchScalarField::mappedFilmPressureFvPatchScalarField
(
    const mappedFilmPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedFilmPressureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Get the mapper and the neighbouring patch
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const fvPatch& patchNbr = mapper.nbrFvPatch();

    // Look up the neighbouring pressure field
    const fvPatchScalarField& pNbr =
        patchNbr.lookupPatchField<volScalarField, scalar>
        (
            internalField().name()
        );

    // Map the neighbouring fluid patch pressure field to this patch
    this->operator==(mapper.fromNeighbour(pNbr));

    // Also assign the mapped pressure to the internal field
    UIndirectList<scalar>
    (
        const_cast<Field<scalar>&>(this->primitiveField()),
        this->patch().faceCells()
    ) = *this;

    zeroGradientFvPatchScalarField::updateCoeffs();
}


void Foam::mappedFilmPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
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
