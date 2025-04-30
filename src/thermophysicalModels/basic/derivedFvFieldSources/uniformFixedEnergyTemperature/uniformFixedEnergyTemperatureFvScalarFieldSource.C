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

#include "uniformFixedEnergyTemperatureFvScalarFieldSource.H"
#include "fvSource.H"
#include "volFields.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformFixedEnergyTemperatureFvScalarFieldSource::
uniformFixedEnergyTemperatureFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    energyCalculatedTemperatureFvScalarFieldSource(iF, dict),
    uniformHe_
    (
        Function1<scalar>::New
        (
            "uniformHe",
            db().time().userUnits(),
            dimEnergy/dimMass,
            dict
        )
    )
{}


Foam::uniformFixedEnergyTemperatureFvScalarFieldSource::
uniformFixedEnergyTemperatureFvScalarFieldSource
(
    const uniformFixedEnergyTemperatureFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    energyCalculatedTemperatureFvScalarFieldSource(field, iF),
    uniformHe_(field.uniformHe_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uniformFixedEnergyTemperatureFvScalarFieldSource::
~uniformFixedEnergyTemperatureFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::uniformFixedEnergyTemperatureFvScalarFieldSource::sourceHeValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return
        DimensionedField<scalar, volMesh>::New
        (
            model.name() + ":" + this->internalField().name() + "SourceHeValue",
            this->internalField().mesh(),
            dimensionedScalar
            (
                this->internalField().dimensions(),
                uniformHe_->value(this->db().time().value())
            )
        );
}


Foam::tmp<Foam::scalarField>
Foam::uniformFixedEnergyTemperatureFvScalarFieldSource::sourceHeValue
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    const scalar v = uniformHe_->value(db().time().value());
    return tmp<scalarField>(new scalarField(source.size(), v));
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::uniformFixedEnergyTemperatureFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return
        DimensionedField<scalar, volMesh>::New
        (
            model.name() + ":" + this->internalField().name() + "InternalCoeff",
            this->internalField().mesh(),
            dimensionedScalar(dimless, scalar(0))
        );
}


Foam::tmp<Foam::scalarField>
Foam::uniformFixedEnergyTemperatureFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return tmp<scalarField>(new scalarField(source.size(), scalar(0)));
}


void Foam::uniformFixedEnergyTemperatureFvScalarFieldSource::write
(
    Ostream& os
) const
{
    fvScalarFieldSource::write(os);
    writeEntry
    (
        os,
        db().time().userUnits(),
        dimEnergy/dimMass,
        uniformHe_()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        uniformFixedEnergyTemperatureFvScalarFieldSource
    );
}

// ************************************************************************* //
