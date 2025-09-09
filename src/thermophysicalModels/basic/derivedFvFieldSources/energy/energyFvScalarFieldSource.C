/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::energyFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const fvScalarFieldSource& Ts = thermo.T().sources()[model.name()];

    if (isA<energyCalculatedTemperatureFvScalarFieldSource>(Ts))
    {
        return
            refCast<const energyCalculatedTemperatureFvScalarFieldSource>(Ts)
           .sourceHeValue(model, source);
    }
    else
    {
        return thermo.he(Ts.sourceValue(model, source), model, source);
    }
}


Foam::tmp<Foam::scalarField>
Foam::energyFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const fvScalarFieldSource& Ts = thermo.T().sources()[model.name()];

    if (isA<energyCalculatedTemperatureFvScalarFieldSource>(Ts))
    {
        return
            refCast<const energyCalculatedTemperatureFvScalarFieldSource>(Ts)
           .sourceHeValue(model, source, cells);
    }
    else
    {
        const scalarField T(Ts.sourceValue(model, source, cells));
        return thermo.he(T, model, source, cells);
    }
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::energyFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const fvScalarFieldSource& Ts = thermo.T().sources()[model.name()];

    return Ts.internalCoeff(model, source);
}


Foam::tmp<Foam::scalarField>
Foam::energyFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const fvScalarFieldSource& Ts = thermo.T().sources()[model.name()];

    return Ts.internalCoeff(model, source, cells);
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
