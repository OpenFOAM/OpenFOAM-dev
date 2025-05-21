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

#include "alphatBoilingWallFunctionFvPatchScalarField.H"
#include "wallBoiling.H"
#include "wallBoilingPhaseChangeRateFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fv::wallBoiling&
Foam::alphatBoilingWallFunctionFvPatchScalarField::model() const
{
    const Foam::fvModels& fvModels =
        Foam::fvModels::New(internalField().mesh());

    if (modelName_ != word::null)
    {
        return refCast<const fv::wallBoiling>(fvModels[modelName_]);
    }

    label wallBoilingFvModeli = -1;

    forAll(fvModels, fvModeli)
    {
        if (!isA<fv::wallBoiling>(fvModels[fvModeli])) continue;

        if (wallBoilingFvModeli != -1)
        {
            FatalErrorInFunction
                << "Multiple wall boiling fvModels found for " << typeName
                << " boundary condition of field " << internalField().name()
                << " on patch " << patch().name() << exit(FatalError);
        }

        wallBoilingFvModeli = fvModeli;
    }

    if (wallBoilingFvModeli == -1)
    {
        FatalErrorInFunction
            << "Wall boiling fvModel not found for " << typeName
            << " boundary condition of field " << internalField().name()
            << " on patch " << patch().name() << exit(FatalError);
    }

    return refCast<const fv::wallBoiling>(fvModels[wallBoilingFvModeli]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * /

Foam::alphatBoilingWallFunctionFvPatchScalarField::
alphatBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    modelName_(dict.lookupOrDefault<word>("model", word::null))
{}


Foam::alphatBoilingWallFunctionFvPatchScalarField::
alphatBoilingWallFunctionFvPatchScalarField
(
    const alphatBoilingWallFunctionFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(psf, p, iF, mapper),
    modelName_(psf.modelName_)
{}


Foam::alphatBoilingWallFunctionFvPatchScalarField::
alphatBoilingWallFunctionFvPatchScalarField
(
    const alphatBoilingWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(psf, iF),
    modelName_(psf.modelName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::alphatBoilingWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fv::wallBoiling& model = this->model();

    const wallBoilingPhaseChangeRateFvPatchScalarField& mDotPf =
        model.mDotPf(patch().index());

    if (&internalField() == &model.alphatLiquid())
    {
        operator==(mDotPf.alphatLiquid());
    }
    else if (&internalField() == &model.alphatVapour())
    {
        operator==(mDotPf.alphatVapour());
    }
    else
    {
        FatalErrorInFunction
            << "Model " << model.name() << " does not provide a turbulent "
            << "thermal diffusivity for phase " << internalField().group()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::alphatBoilingWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);

    writeEntryIfDifferent(os, "model", word::null, modelName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        alphatBoilingWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
