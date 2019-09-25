/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    const fvPatchFieldMapper& mapper
)
:
    calculatedFvPatchScalarField(ptf, p, iF, mapper),
    heGradient_(mapper(ptf.heGradient_))
{}


Foam::gradientEnergyCalculatedTemperatureFvPatchScalarField::
gradientEnergyCalculatedTemperatureFvPatchScalarField
(
    const gradientEnergyCalculatedTemperatureFvPatchScalarField& ptf
)
:
    calculatedFvPatchScalarField(ptf),
    heGradient_(ptf.heGradient_)
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

void Foam::gradientEnergyCalculatedTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    calculatedFvPatchScalarField::autoMap(m);
    m(heGradient_, heGradient_);
}


void Foam::gradientEnergyCalculatedTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    calculatedFvPatchScalarField::rmap(ptf, addr);

    const gradientEnergyCalculatedTemperatureFvPatchScalarField& mptf =
        refCast<const gradientEnergyCalculatedTemperatureFvPatchScalarField>
        (
            ptf
        );

    heGradient_.rmap(mptf.heGradient_, addr);
}


// ************************************************************************* //
