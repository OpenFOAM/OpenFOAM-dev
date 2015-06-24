/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "gradientUnburntEnthalpyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "psiuReactionThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradientUnburntEnthalpyFvPatchScalarField::
gradientUnburntEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


Foam::gradientUnburntEnthalpyFvPatchScalarField::
gradientUnburntEnthalpyFvPatchScalarField
(
    const gradientUnburntEnthalpyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::gradientUnburntEnthalpyFvPatchScalarField::
gradientUnburntEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{}


Foam::gradientUnburntEnthalpyFvPatchScalarField::
gradientUnburntEnthalpyFvPatchScalarField
(
    const gradientUnburntEnthalpyFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf)
{}


Foam::gradientUnburntEnthalpyFvPatchScalarField::
gradientUnburntEnthalpyFvPatchScalarField
(
    const gradientUnburntEnthalpyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradientUnburntEnthalpyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const psiuReactionThermo& thermo = db().lookupObject<psiuReactionThermo>
    (
        basicThermo::dictName
    );

    const label patchi = patch().index();

    const scalarField& pw = thermo.p().boundaryField()[patchi];
    fvPatchScalarField& Tw =
        const_cast<fvPatchScalarField&>(thermo.Tu().boundaryField()[patchi]);

    Tw.evaluate();

    gradient() = thermo.Cp(pw, Tw, patchi)*Tw.snGrad()
      + patch().deltaCoeffs()*
        (
            thermo.heu(pw, Tw, patchi)
          - thermo.heu(pw, Tw, patch().faceCells())
        );

    fixedGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        gradientUnburntEnthalpyFvPatchScalarField
    );
}

// ************************************************************************* //
