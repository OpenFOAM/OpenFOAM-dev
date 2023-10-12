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

#include "uniformInletOutletEnergyTemperatureFvScalarFieldSource.H"
#include "fvSource.H"
#include "volFields.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformInletOutletEnergyTemperatureFvScalarFieldSource::
uniformInletOutletEnergyTemperatureFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    energyCalculatedTemperatureFvScalarFieldSource(iF, dict),
    uniformInletHe_(Function1<scalar>::New("uniformInletHe", dict))
{}


Foam::uniformInletOutletEnergyTemperatureFvScalarFieldSource::
uniformInletOutletEnergyTemperatureFvScalarFieldSource
(
    const uniformInletOutletEnergyTemperatureFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    energyCalculatedTemperatureFvScalarFieldSource(field, iF),
    uniformInletHe_(field.uniformInletHe_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uniformInletOutletEnergyTemperatureFvScalarFieldSource::
~uniformInletOutletEnergyTemperatureFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::uniformInletOutletEnergyTemperatureFvScalarFieldSource::sourceHeValue
(
    const fvSource& source
) const
{
    const scalar t = this->db().time().userTimeValue();
    const scalar v = uniformInletHe_->value(t);
    return tmp<scalarField>(new scalarField(source.nCells(), v));
}


Foam::tmp<Foam::scalarField>
Foam::uniformInletOutletEnergyTemperatureFvScalarFieldSource::internalCoeff
(
    const fvSource& source
) const
{
    return
        neg0(source.source(internalField().name()))
       *scalarField(source.nCells(), scalar(1));
}


void Foam::uniformInletOutletEnergyTemperatureFvScalarFieldSource::write
(
    Ostream& os
) const
{
    fvScalarFieldSource::write(os);
    writeEntry(os, uniformInletHe_());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        uniformInletOutletEnergyTemperatureFvScalarFieldSource
    );
}

// ************************************************************************* //
