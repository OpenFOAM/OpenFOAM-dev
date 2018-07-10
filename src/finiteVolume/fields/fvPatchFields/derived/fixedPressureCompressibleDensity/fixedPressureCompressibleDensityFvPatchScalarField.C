/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "fixedPressureCompressibleDensityFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedPressureCompressibleDensityFvPatchScalarField::
fixedPressureCompressibleDensityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    pName_("p")
{}


Foam::fixedPressureCompressibleDensityFvPatchScalarField::
fixedPressureCompressibleDensityFvPatchScalarField
(
    const fixedPressureCompressibleDensityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    pName_(ptf.pName_)
{}


Foam::fixedPressureCompressibleDensityFvPatchScalarField::
fixedPressureCompressibleDensityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    pName_(dict.lookupOrDefault<word>("p", "p"))
{}


Foam::fixedPressureCompressibleDensityFvPatchScalarField::
fixedPressureCompressibleDensityFvPatchScalarField
(
    const fixedPressureCompressibleDensityFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    pName_(ptf.pName_)
{}


Foam::fixedPressureCompressibleDensityFvPatchScalarField::
fixedPressureCompressibleDensityFvPatchScalarField
(
    const fixedPressureCompressibleDensityFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    pName_(ptf.pName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedPressureCompressibleDensityFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& pp =
        patch().lookupPatchField<volScalarField, scalar>(pName_);

    const dictionary& thermoProps =
        db().lookupObject<IOdictionary>("thermodynamicProperties");

    const scalar rholSat =
        dimensionedScalar(thermoProps.lookup("rholSat")).value();

    const scalar pSat =
        dimensionedScalar(thermoProps.lookup("pSat")).value();

    const scalar psil = dimensionedScalar(thermoProps.lookup("psil")).value();

    operator==(rholSat + psil*(pp - pSat));

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::fixedPressureCompressibleDensityFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedPressureCompressibleDensityFvPatchScalarField
    );
}

// ************************************************************************* //
