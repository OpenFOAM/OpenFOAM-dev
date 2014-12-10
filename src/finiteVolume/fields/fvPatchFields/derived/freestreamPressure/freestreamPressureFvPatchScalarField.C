/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "freestreamPressureFvPatchScalarField.H"
#include "freestreamFvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF)
{}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict)
{}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& wbppsf
)
:
    zeroGradientFvPatchScalarField(wbppsf)
{}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freestreamPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const freestreamFvPatchVectorField& Up =
        refCast<const freestreamFvPatchVectorField>
        (
            patch().lookupPatchField<volVectorField, vector>("U")
        );

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>("phi");

    fvsPatchField<scalar>& phip =
        const_cast<fvsPatchField<scalar>&>
        (
            patch().patchField<surfaceScalarField, scalar>(phi)
        );

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        phip = patch().Sf() & Up.freestreamValue();
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>("rho");

        phip = rhop*(patch().Sf() & Up.freestreamValue());
    }
    else
    {
        FatalErrorIn("freestreamPressureFvPatchScalarField::updateCoeffs()")
            << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    zeroGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        freestreamPressureFvPatchScalarField
    );
}

// ************************************************************************* //
