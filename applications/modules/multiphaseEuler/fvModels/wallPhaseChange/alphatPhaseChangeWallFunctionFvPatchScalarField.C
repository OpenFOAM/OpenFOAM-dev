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

#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"
#include "wallPhaseChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::UPtrList<const Foam::fv::wallPhaseChange>&
Foam::alphatPhaseChangeWallFunctionFvPatchScalarField::models() const
{
    if (models_.size()) return models_;

    const Foam::fvModels& fvModels =
        Foam::fvModels::New(internalField().mesh());

    forAll(fvModels, fvModeli)
    {
        if (!isA<fv::wallPhaseChange>(fvModels[fvModeli])) continue;

        const fv::wallPhaseChange& model =
            refCast<const fv::wallPhaseChange>(fvModels[fvModeli]);

        if
        (
            &internalField() == &model.alphats().first()
         || &internalField() == &model.alphats().second()
        )
        {
            models_.append(&model);
        }
    }

    if (models_.empty())
    {
        FatalErrorInFunction
            << "No wall phase-change models found for " << typeName
            << " boundary condition of field " << internalField().name()
            << " on patch " << patch().name() << exit(FatalError);
    }

    return models_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * /

Foam::alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    models_()
{}


Foam::alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const alphatPhaseChangeWallFunctionFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(psf, p, iF, mapper),
    models_()
{}


Foam::alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const alphatPhaseChangeWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(psf, iF),
    models_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::alphatPhaseChangeWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const UPtrList<const fv::wallPhaseChange>& models = this->models();

    // Apply the first model to the entire patch, and subsequent models only to
    // the active portion of the patch. This assumes that the models active
    // regions do not overlap. If it turns out that we do need to support
    // multiple models simultaneously causing phase-change on the same face,
    // then these models will need to get more collaborative and we have to
    // think about how to do that.
    forAll(models, modeli)
    {
        const fv::wallPhaseChange& model = models[modeli];

        const scalarField& alphat =
            &internalField() == &model.alphats().first()
          ? model.alphats(patch().index()).first()
          : model.alphats(patch().index()).second();

        if (modeli == 0)
        {
            operator==(alphat);
        }
        else
        {
            const scalarField& active = model.active(patch().index());
            operator==((1 - active)*(*this) + active*alphat);
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        alphatPhaseChangeWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
