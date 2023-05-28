/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "uniformDensityHydrostaticPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    rho_(dict.lookup<scalar>("rhoRef")),
    pRef_(dict.lookup<scalar>("pRef")),
    pRefPointSpecified_(dict.found("pRefPoint")),
    pRefPoint_(dict.lookupOrDefault<vector>("pRefPoint", Zero))
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


Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
    const uniformDensityHydrostaticPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    rho_(ptf.rho_),
    pRef_(ptf.pRef_),
    pRefPointSpecified_(ptf.pRefPointSpecified_),
    pRefPoint_(ptf.pRefPoint_)
{}


Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
    const uniformDensityHydrostaticPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    rho_(ptf.rho_),
    pRef_(ptf.pRef_),
    pRefPointSpecified_(ptf.pRefPointSpecified_),
    pRefPoint_(ptf.pRefPoint_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniformDensityHydrostaticPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    scalar ghRef = g.value() & pRefPoint_;

    if (!pRefPointSpecified_)
    {
        const uniformDimensionedScalarField& hRef =
            db().lookupObject<uniformDimensionedScalarField>("hRef");

        ghRef = - mag(g.value())*hRef.value();
    }

    operator==
    (
        pRef_
      + rho_*((g.value() & patch().Cf()) - ghRef)
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::uniformDensityHydrostaticPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "rhoRef", rho_);
    writeEntry(os, "pRef", pRef_);
    if (pRefPointSpecified_)
    {
        writeEntry(os, "pRefPoint", pRefPoint_);
    }
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        uniformDensityHydrostaticPressureFvPatchScalarField
    );
}

// ************************************************************************* //
