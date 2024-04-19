/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "prghCyclicPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prghCyclicPressureFvPatchScalarField::prghCyclicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    jumpCyclicFvPatchScalarField(p, iF, dict),
    rhoName_
    (
        cyclicPatch().owner()
      ? dict.lookupOrDefault<word>("rho", "rho")
      : word::null
    ),
    rhoInf_
    (
        cyclicPatch().owner()
      ? dict.lookup<scalar>("rhoInf", dimDensity)
      : NaN
    ),
    jump_(p.size(), Zero)
{
    if (dict.found("jump"))
    {
        jump_ = scalarField("jump", iF.dimensions(), dict, p.size());
    }

    evaluateNoUpdateCoeffs();
}


Foam::prghCyclicPressureFvPatchScalarField::prghCyclicPressureFvPatchScalarField
(
    const prghCyclicPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    jumpCyclicFvPatchScalarField(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_),
    rhoInf_(ptf.rhoInf_),
    jump_(mapper(ptf.jump_))
{}


Foam::prghCyclicPressureFvPatchScalarField::prghCyclicPressureFvPatchScalarField
(
    const prghCyclicPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    jumpCyclicFvPatchScalarField(ptf, iF),
    rhoName_(ptf.rhoName_),
    rhoInf_(ptf.rhoInf_),
    jump_(ptf.jump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::prghCyclicPressureFvPatchScalarField::jump() const
{
    return jump_;
}


void Foam::prghCyclicPressureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    jumpCyclicFvPatchScalarField::map(ptf, mapper);

    const prghCyclicPressureFvPatchScalarField& tiptf =
        refCast<const prghCyclicPressureFvPatchScalarField>(ptf);

    mapper(jump_, tiptf.jump_);
}


void Foam::prghCyclicPressureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    jumpCyclicFvPatchScalarField::reset(ptf);

    const prghCyclicPressureFvPatchScalarField& tiptf =
        refCast<const prghCyclicPressureFvPatchScalarField>(ptf);

    jump_.reset(tiptf.jump_);
}


void Foam::prghCyclicPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const volScalarField& vf =
        static_cast<const volScalarField&>(internalField());

    const prghCyclicPressureFvPatchScalarField& nbrPf =
        refCast<const prghCyclicPressureFvPatchScalarField>(nbrPatchField());

    const label patchi = patch().index(), nbrPatchi = nbrPf.patch().index();

    // Buoyancy fields
    const volScalarField& rhoVf =
        db().lookupObject<volScalarField>
        (cyclicPatch().owner() ? rhoName_ : nbrPf.rhoName_);
    const volScalarField::Boundary& rhoBf = rhoVf.boundaryField();
    const surfaceScalarField::Boundary& ghfBf =
        db().lookupObject<surfaceScalarField>("ghf").boundaryField();

    // Pressure solution fields
    const surfaceScalarField::Boundary& rAUfBf =
        db().lookupObject<surfaceScalarField>("rAUf").boundaryField();
    const surfaceScalarField::Boundary& phiHbyABf =
        db().lookupObject<surfaceScalarField>("phiHbyA").boundaryField();

    // Delta coefficients for this field
    const tmp<surfaceScalarField> deltaCoeffsSf =
        fv::snGradScheme<scalar>::New
        (
            vf.mesh(),
            vf.mesh().schemes().snGrad(internalField().name())
        )->deltaCoeffs(vf);
    const scalarField& deltaCoeffsPf =
        deltaCoeffsSf->boundaryField()[patch().index()];

    // Calculate the jump
    jump_ =
        (rhoBf[patchi] - (cyclicPatch().owner() ? rhoInf_ : nbrPf.rhoInf_))
       *(ghfBf[nbrPatchi] - ghfBf[patchi])
      + (
            phiHbyABf[patchi]/rAUfBf[patchi]/patch().magSf()
          + phiHbyABf[nbrPatchi]/rAUfBf[nbrPatchi]/nbrPf.patch().magSf()
        )*patch().weights()/deltaCoeffsPf;

    jumpCyclicFvPatchScalarField::updateCoeffs();
}


void Foam::prghCyclicPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    if (cyclicPatch().owner())
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
        writeEntry(os, "rhoInf", rhoInf_);
    }

    writeEntry(os, "jump", jump_);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        prghCyclicPressureFvPatchScalarField
    );
}

// ************************************************************************* //
