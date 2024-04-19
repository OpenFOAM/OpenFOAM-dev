/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "dimensionSet.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

HashTable<dimensionSet>* dimensionsPtr_(nullptr);

// Delete the above data at the end of the run
struct deleteDimensionsPtr
{
    ~deleteDimensionsPtr()
    {
        deleteDemandDrivenData(dimensionsPtr_);
    }
};

deleteDimensionsPtr deleteDimensionsPtr_;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::dimensionSet Foam::dimless(0, 0, 0, 0, 0, 0, 0);

const Foam::dimensionSet Foam::dimMass(1, 0, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimLength(0, 1, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimTime(0, 0, 1, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimTemperature(0, 0, 0, 1, 0, 0, 0);
const Foam::dimensionSet Foam::dimMoles(0, 0, 0, 0, 1, 0, 0);
const Foam::dimensionSet Foam::dimCurrent(0, 0, 0, 0, 0, 1, 0);
const Foam::dimensionSet Foam::dimLuminousIntensity(0, 0, 0, 0, 0, 0, 1);

const Foam::dimensionSet Foam::dimArea(sqr(dimLength));
const Foam::dimensionSet Foam::dimVolume(pow3(dimLength));

const Foam::dimensionSet Foam::dimRate(dimless/dimTime);

const Foam::dimensionSet Foam::dimVelocity(dimLength/dimTime);
const Foam::dimensionSet Foam::dimMomentum(dimMass*dimVelocity);
const Foam::dimensionSet Foam::dimAcceleration(dimVelocity/dimTime);

const Foam::dimensionSet Foam::dimDensity(dimMass/dimVolume);
const Foam::dimensionSet Foam::dimForce(dimMass*dimAcceleration);
const Foam::dimensionSet Foam::dimEnergy(dimForce*dimLength);
const Foam::dimensionSet Foam::dimPower(dimEnergy/dimTime);

const Foam::dimensionSet Foam::dimPressure(dimForce/dimArea);
const Foam::dimensionSet Foam::dimKinematicPressure(dimPressure/dimDensity);
const Foam::dimensionSet Foam::dimCompressibility(dimDensity/dimPressure);
const Foam::dimensionSet Foam::dimGasConstant(dimEnergy/dimMass/dimTemperature);
const Foam::dimensionSet Foam::dimSpecificHeatCapacity(dimGasConstant);
const Foam::dimensionSet Foam::dimKinematicViscosity(dimArea/dimTime);
const Foam::dimensionSet Foam::dimDynamicViscosity
(
    dimDensity*dimKinematicViscosity
);
const Foam::dimensionSet Foam::dimThermalConductivity
(
    dimPower/dimLength/dimTemperature
);

const Foam::dimensionSet Foam::dimVolumetricFlux(dimArea*dimVelocity);
const Foam::dimensionSet Foam::dimMassFlux(dimDensity*dimVolumetricFlux);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::HashTable<Foam::dimensionSet>& Foam::dimensions()
{
    static const NamedEnum<dimensionSet::dimensionType, 7>& names =
        dimensionSet::dimensionTypeNames_;

    if (!dimensionsPtr_)
    {
        dimensionsPtr_ = new HashTable<dimensionSet>();

        dimensionsPtr_->insert(names[dimensionSet::MASS], dimMass);
        dimensionsPtr_->insert(names[dimensionSet::LENGTH], dimLength);
        dimensionsPtr_->insert(names[dimensionSet::TIME], dimTime);
        dimensionsPtr_->insert
        (
            names[dimensionSet::TEMPERATURE],
            dimTemperature
        );
        dimensionsPtr_->insert(names[dimensionSet::MOLES], dimMoles);
        dimensionsPtr_->insert(names[dimensionSet::CURRENT], dimCurrent);
        dimensionsPtr_->insert
        (
            names[dimensionSet::LUMINOUS_INTENSITY],
            dimLuminousIntensity
        );

        dimensionsPtr_->insert("area", dimArea);
        dimensionsPtr_->insert("volume", dimVolume);

        dimensionsPtr_->insert("rate", dimRate);

        dimensionsPtr_->insert("velocity", dimVelocity);
        dimensionsPtr_->insert("momentum", dimMomentum);
        dimensionsPtr_->insert("acceleration", dimAcceleration);

        dimensionsPtr_->insert("density", dimDensity);
        dimensionsPtr_->insert("force", dimForce);
        dimensionsPtr_->insert("energy", dimEnergy);
        dimensionsPtr_->insert("power", dimPower);

        dimensionsPtr_->insert("pressure", dimPressure);
        dimensionsPtr_->insert("kinematicPressure", dimKinematicPressure);
        dimensionsPtr_->insert("compressibility", dimCompressibility);
        dimensionsPtr_->insert("gasConstant", dimGasConstant);
        dimensionsPtr_->insert("specificHeatCapacity", dimSpecificHeatCapacity);
        dimensionsPtr_->insert("kinematicViscosity", dimKinematicViscosity);
        dimensionsPtr_->insert("dynamicViscosity", dimDynamicViscosity);
        dimensionsPtr_->insert("thermalConductivity", dimThermalConductivity);

        dimensionsPtr_->insert("volumetricFlux", dimVolumetricFlux);
        dimensionsPtr_->insert("massFlux", dimMassFlux);
    }

    return *dimensionsPtr_;
}


// ************************************************************************* //
