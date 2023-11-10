/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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

#include "gradientEnergyCalculatedTemperatureFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        gradientEnergyCalculatedTemperatureFvPatchScalarField,
        0
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradientEnergyCalculatedTemperatureFvPatchScalarField::
gradientEnergyCalculatedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(p, iF),
    heGradient_(p.size())
{}


Foam::gradientEnergyCalculatedTemperatureFvPatchScalarField::
gradientEnergyCalculatedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    calculatedFvPatchScalarField(p, iF),
    heGradient_(p.size())
{
    calculatedFvPatchScalarField::evaluate();
}


Foam::gradientEnergyCalculatedTemperatureFvPatchScalarField::
gradientEnergyCalculatedTemperatureFvPatchScalarField
(
    const gradientEnergyCalculatedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    calculatedFvPatchScalarField(ptf, p, iF, mapper),
    heGradient_(mapper(ptf.heGradient_))
{}


Foam::gradientEnergyCalculatedTemperatureFvPatchScalarField::
gradientEnergyCalculatedTemperatureFvPatchScalarField
(
    const gradientEnergyCalculatedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(ptf, iF),
    heGradient_(ptf.heGradient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradientEnergyCalculatedTemperatureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    calculatedFvPatchScalarField::map(ptf, mapper);

    const gradientEnergyCalculatedTemperatureFvPatchScalarField& mptf =
        refCast<const gradientEnergyCalculatedTemperatureFvPatchScalarField>
        (
            ptf
        );

    mapper(heGradient_, mptf.heGradient_);
}


void Foam::gradientEnergyCalculatedTemperatureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    calculatedFvPatchScalarField::reset(ptf);

    const gradientEnergyCalculatedTemperatureFvPatchScalarField& mptf =
        refCast<const gradientEnergyCalculatedTemperatureFvPatchScalarField>
        (
            ptf
        );

    heGradient_.reset(mptf.heGradient_);
}


// ************************************************************************* //
