/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "energyFvScalarFieldSource.H"
#include "energyCalculatedTemperatureFvScalarFieldSource.H"
#include "fvSource.H"
#include "volFields.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyFvScalarFieldSource::energyFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(iF)
{}


Foam::energyFvScalarFieldSource::energyFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict)
{}


Foam::energyFvScalarFieldSource::energyFvScalarFieldSource
(
    const energyFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::energyFvScalarFieldSource::~energyFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::energyFvScalarFieldSource::sourceValue(const fvSource& source) const
{
    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const fvScalarFieldSource& Ts = thermo.T().sources()[source.name()];

    if (isA<energyCalculatedTemperatureFvScalarFieldSource>(Ts))
    {
        return
            refCast<const energyCalculatedTemperatureFvScalarFieldSource>(Ts)
           .sourceHeValue(source);
    }
    else
    {
        return thermo.he(Ts.sourceValue(source), source);
    }
}


Foam::tmp<Foam::scalarField>
Foam::energyFvScalarFieldSource::internalCoeff(const fvSource& source) const
{
    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const fvScalarFieldSource& Ts = thermo.T().sources()[source.name()];

    return Ts.internalCoeff(source);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeNullConstructableTypeFieldSource
    (
        fvScalarFieldSource,
        energyFvScalarFieldSource
    );
}

// ************************************************************************* //
