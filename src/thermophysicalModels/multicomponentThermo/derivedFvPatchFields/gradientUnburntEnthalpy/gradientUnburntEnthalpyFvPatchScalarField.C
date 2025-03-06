/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
#include "fieldMapper.H"
#include "psiuMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradientUnburntEnthalpyFvPatchScalarField::
gradientUnburntEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    gradient() = scalar(0);
}


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
    const gradientUnburntEnthalpyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper, false)
{
    map(ptf, mapper);
}


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

void Foam::gradientUnburntEnthalpyFvPatchScalarField::map
(
    const gradientUnburntEnthalpyFvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    // Unmapped faces are considered zero-gradient/adiabatic
    // until they are corrected later
    mapper(*this, ptf, [&](){ return this->patchInternalField(); });
    mapper(gradient(), ptf.gradient(), scalar(0));
}


void Foam::gradientUnburntEnthalpyFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    map(refCast<const gradientUnburntEnthalpyFvPatchScalarField>(ptf), mapper);
}


void Foam::gradientUnburntEnthalpyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const psiuMulticomponentThermo& thermo =
        db().lookupObject<psiuMulticomponentThermo>
        (
            physicalProperties::typeName
        );

    const label patchi = patch().index();

    fvPatchScalarField& Tw =
        const_cast<fvPatchScalarField&>(thermo.Tu().boundaryField()[patchi]);

    Tw.evaluate();

    gradient() = thermo.Cp(Tw, patchi)*Tw.snGrad()
      + patch().deltaCoeffs()*
        (
            thermo.heu(Tw, patchi)
          - thermo.heu(Tw, patch().faceCells())
        );

    fixedGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeNullConstructablePatchTypeField
    (
        fvPatchScalarField,
        gradientUnburntEnthalpyFvPatchScalarField
    );
}

// ************************************************************************* //
