/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "prghUniformDensityHydrostaticPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prghUniformDensityHydrostaticPressureFvPatchScalarField::
prghUniformDensityHydrostaticPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    pRef_(0),
    rhoRef_(0),
    rhoName_("rho")
{}


Foam::prghUniformDensityHydrostaticPressureFvPatchScalarField::
prghUniformDensityHydrostaticPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    pRef_(readScalar(dict.lookup("pRef"))),
    rhoRef_(readScalar(dict.lookup("rhoRef"))),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(pRef_);
    }
}


Foam::prghUniformDensityHydrostaticPressureFvPatchScalarField::
prghUniformDensityHydrostaticPressureFvPatchScalarField
(
    const prghUniformDensityHydrostaticPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    pRef_(ptf.pRef_),
    rhoRef_(ptf.rhoRef_),
    rhoName_(ptf.rhoName_)
{}


Foam::prghUniformDensityHydrostaticPressureFvPatchScalarField::
prghUniformDensityHydrostaticPressureFvPatchScalarField
(
    const prghUniformDensityHydrostaticPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    pRef_(ptf.pRef_),
    rhoRef_(ptf.rhoRef_),
    rhoName_(ptf.rhoName_)
{}


Foam::prghUniformDensityHydrostaticPressureFvPatchScalarField::
prghUniformDensityHydrostaticPressureFvPatchScalarField
(
    const prghUniformDensityHydrostaticPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    pRef_(ptf.pRef_),
    rhoRef_(ptf.rhoRef_),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::prghUniformDensityHydrostaticPressureFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& rhop = patch().lookupPatchField<volScalarField, scalar>
    (
        rhoName_
    );

    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    const uniformDimensionedScalarField& hRef =
        db().lookupObject<uniformDimensionedScalarField>("hRef");

    dimensionedScalar ghRef
    (
        mag(g.value()) > small
      ? g & (cmptMag(g.value())/mag(g.value()))*hRef
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    operator==
    (
        pRef_ - (rhop - rhoRef_)*((g.value() & patch().Cf()) - ghRef.value())
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::prghUniformDensityHydrostaticPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("pRef") << pRef_ << token::END_STATEMENT << nl;
    os.writeKeyword("rhoRef") << rhoRef_ << token::END_STATEMENT << nl;
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        prghUniformDensityHydrostaticPressureFvPatchScalarField
    );
}

// ************************************************************************* //
