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

#include "alphatCondensationWallFunctionFvPatchScalarField.H"
#include "wallCondensation.H"
#include "wallCondensationPhaseChangeRateFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fv::wallCondensation&
Foam::alphatCondensationWallFunctionFvPatchScalarField::model() const
{
    const Foam::fvModels& fvModels =
        Foam::fvModels::New(internalField().mesh());

    if (modelName_ != word::null)
    {
        return refCast<const fv::wallCondensation>(fvModels[modelName_]);
    }

    label wallCondensationFvModeli = -1;

    forAll(fvModels, fvModeli)
    {
        if (!isA<fv::wallCondensation>(fvModels[fvModeli])) continue;

        if (wallCondensationFvModeli != -1)
        {
            FatalErrorInFunction
                << "Multiple wall condensation fvModels found for " << typeName
                << " boundary condition of field " << internalField().name()
                << " on patch " << patch().name() << exit(FatalError);
        }

        wallCondensationFvModeli = fvModeli;
    }

    if (wallCondensationFvModeli == -1)
    {
        FatalErrorInFunction
            << "Wall condensation fvModel not found for " << typeName
            << " boundary condition of field " << internalField().name()
            << " on patch " << patch().name() << exit(FatalError);
    }

    return
        refCast<const fv::wallCondensation>(fvModels[wallCondensationFvModeli]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * /

Foam::alphatCondensationWallFunctionFvPatchScalarField::
alphatCondensationWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    modelName_(dict.lookupOrDefault<word>("model", word::null))
{}


Foam::alphatCondensationWallFunctionFvPatchScalarField::
alphatCondensationWallFunctionFvPatchScalarField
(
    const alphatCondensationWallFunctionFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(psf, p, iF, mapper),
    modelName_(psf.modelName_)
{}


Foam::alphatCondensationWallFunctionFvPatchScalarField::
alphatCondensationWallFunctionFvPatchScalarField
(
    const alphatCondensationWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(psf, iF),
    modelName_(psf.modelName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::alphatCondensationWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fv::wallCondensation& model = this->model();

    const wallCondensationPhaseChangeRateFvPatchScalarField& mDotPf =
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


void
Foam::alphatCondensationWallFunctionFvPatchScalarField::write(Ostream& os) const
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
        alphatCondensationWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
