/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::dimensionSet Foam::dimensions::invalid
(
    vGreat, vGreat, vGreat, vGreat, vGreat, vGreat, vGreat
);

const Foam::dimensionSet Foam::dimensions::dimless(0, 0, 0, 0, 0, 0, 0);

Foam::HashTable<Foam::dimensionSet> Foam::dimensions::table;

const Foam::dimensionSet Foam::dimensions::mass
(
    Foam::dimensionSet::dimensionTypeNames_[Foam::dimensionSet::MASS],
    1, 0, 0, 0, 0, 0, 0
);

const Foam::dimensionSet Foam::dimensions::length
(
    Foam::dimensionSet::dimensionTypeNames_[Foam::dimensionSet::LENGTH],
    0, 1, 0, 0, 0, 0, 0
);

const Foam::dimensionSet Foam::dimensions::time
(
    Foam::dimensionSet::dimensionTypeNames_[Foam::dimensionSet::TIME],
    0, 0, 1, 0, 0, 0, 0
);

const Foam::dimensionSet Foam::dimensions::temperature
(
    Foam::dimensionSet::dimensionTypeNames_ [Foam::dimensionSet::TEMPERATURE],
    0, 0, 0, 1, 0, 0, 0
);

const Foam::dimensionSet Foam::dimensions::moles
(
    Foam::dimensionSet::dimensionTypeNames_[Foam::dimensionSet::MOLES],
    0, 0, 0, 0, 1, 0, 0
);

const Foam::dimensionSet Foam::dimensions::current
(
    Foam::dimensionSet::dimensionTypeNames_[Foam::dimensionSet::CURRENT],
    0, 0, 0, 0, 0, 1, 0
);

const Foam::dimensionSet Foam::dimensions::luminousIntensity
(
    Foam::dimensionSet::dimensionTypeNames_
    [
        Foam::dimensionSet::LUMINOUS_INTENSITY
    ],
    0, 0, 0, 0, 0, 0, 1
);

const Foam::dimensionSet Foam::dimensions::area("area", sqr(length));
const Foam::dimensionSet Foam::dimensions::volume("volume", pow3(length));

const Foam::dimensionSet Foam::dimensions::rate("rate", dimless/time);

const Foam::dimensionSet Foam::dimensions::velocity("velocity", length/time);
const Foam::dimensionSet Foam::dimensions::momentum("momentum", mass*velocity);
const Foam::dimensionSet Foam::dimensions::acceleration
(
    "acceleration",
    velocity/time
);
const Foam::dimensionSet Foam::dimensions::force("force", mass*acceleration);

const Foam::dimensionSet Foam::dimensions::density("density", mass/volume);
const Foam::dimensionSet Foam::dimensions::momentumDensity
(
    "momentumDensity",
    momentum/volume
);
const Foam::dimensionSet Foam::dimensions::forceDensity
(
    "forceDensity",
    force/volume
);

const Foam::dimensionSet Foam::dimensions::energy("energy", force*length);
const Foam::dimensionSet Foam::dimensions::energyDensity
(
    "energyDensity",
    energy/volume
);
const Foam::dimensionSet Foam::dimensions::specificEnergy
(
    "specificEnergy",
    energy/mass
);

const Foam::dimensionSet Foam::dimensions::kineticEnergy
(
    "kineticEnergy",
    energy
);
const Foam::dimensionSet Foam::dimensions::kineticEnergyDensity
(
    "kineticEnergyDensity",
    kineticEnergy/volume
);

const Foam::dimensionSet Foam::dimensions::power("power", energy/time);
const Foam::dimensionSet Foam::dimensions::powerDensity
(
    "powerDensity",
    power/volume
);
const Foam::dimensionSet Foam::dimensions::specificPower
(
    "specificPower",
    power/mass
);

const Foam::dimensionSet Foam::dimensions::entropy
(
    "entropy",
    energy/temperature
);
const Foam::dimensionSet Foam::dimensions::specificEntropy
(
    "specificEntropy",
    entropy/mass
);

const Foam::dimensionSet Foam::dimensions::heatCapacity
(
    "heatCapacity",
    energy/temperature
);
const Foam::dimensionSet Foam::dimensions::specificHeatCapacity
(
    "specificHeatCapacity",
    heatCapacity/mass
);

const Foam::dimensionSet Foam::dimensions::gasConstant
(
    "gasConstant",
    specificHeatCapacity
);

const Foam::dimensionSet Foam::dimensions::pressure("pressure", force/area);
const Foam::dimensionSet Foam::dimensions::kinematicPressure
(
    "kinematicPressure",
    pressure/density
);
const Foam::dimensionSet Foam::dimensions::compressibility
(
    "compressibility",
    density/pressure
);

const Foam::dimensionSet Foam::dimensions::kinematicViscosity
(
    "kinematicViscosity",
    area/time
);
const Foam::dimensionSet Foam::dimensions::dynamicViscosity
(
    "dynamicViscosity",
    density*kinematicViscosity
);
const Foam::dimensionSet Foam::dimensions::kinematicDiffusivity
(
    "kinematicDiffusivity",
    kinematicViscosity
);
const Foam::dimensionSet Foam::dimensions::dynamicDiffusivity
(
    "dynamicDiffusivity",
    dynamicViscosity
);
const Foam::dimensionSet Foam::dimensions::thermalConductivity
(
    "thermalConductivity",
    power/length/temperature
);

const Foam::dimensionSet Foam::dimensions::turbulentKineticEnergy
(
    "turbulentKineticEnergy",
    sqr(velocity)
);
const Foam::dimensionSet Foam::dimensions::kinematicStress
(
    "kinematicStress",
    turbulentKineticEnergy
);
const Foam::dimensionSet Foam::dimensions::ReynoldsStress
(
    "ReynoldsStress",
    kinematicStress
);
const Foam::dimensionSet Foam::dimensions::turbulentEpsilon
(
    "turbulentEpsilon",
    turbulentKineticEnergy/time
);
const Foam::dimensionSet Foam::dimensions::turbulentOmega
(
    "turbulentOmega",
    rate
);
const Foam::dimensionSet Foam::dimensions::turbulentViscosity
(
    "turbulentViscosity",
    kinematicViscosity
);

const Foam::dimensionSet Foam::dimensions::volumetricFlux
(
    "volumetricFlux",
    area*velocity
);
const Foam::dimensionSet Foam::dimensions::volumetricFluxDensity
(
    "volumetricFluxDensity",
    volumetricFlux/area
);

const Foam::dimensionSet Foam::dimensions::massFlux
(
    "massFlux",
    density*volumetricFlux
);
const Foam::dimensionSet Foam::dimensions::massFluxDensity
(
    "massFluxDensity",
    massFlux/area
);

const Foam::dimensionSet Foam::dimensions::heatFlux("heatFlux", power);
const Foam::dimensionSet Foam::dimensions::heatFluxDensity
(
    "heatFluxDensity",
    heatFlux/area
);

const Foam::dimensionSet Foam::dimensions::charge
(
    "charge",
    current*time
);
const Foam::dimensionSet Foam::dimensions::chargeDensity
(
    "chargeDensity",
    charge/volume
);
const Foam::dimensionSet Foam::dimensions::electricPotential
(
    "electricPotential",
    power/current
);
const Foam::dimensionSet Foam::dimensions::electricalConductivity
(
    "electricalConductivity",
    sqr(dimensions::current)/dimensions::length/dimensions::power
);

const Foam::dimensionSet Foam::dimensions::magneticFluxDensity
(
    "magneticFluxDensity",
    force/(length*current)
);
const Foam::dimensionSet Foam::dimensions::magneticFluxPressure
(
    "magneticFluxPressure",
    magneticFluxDensity*velocity
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::dimensionSet& Foam::dimMass = dimensions::mass;
const Foam::dimensionSet& Foam::dimLength = dimensions::length;
const Foam::dimensionSet& Foam::dimTime = dimensions::time;
const Foam::dimensionSet& Foam::dimTemperature = dimensions::temperature;
const Foam::dimensionSet& Foam::dimMoles = dimensions::moles;

const Foam::dimensionSet& Foam::dimArea = dimensions::area;
const Foam::dimensionSet& Foam::dimVolume = dimensions::volume;

const Foam::dimensionSet& Foam::dimRate = dimensions::rate;

const Foam::dimensionSet& Foam::dimVelocity = dimensions::velocity;
const Foam::dimensionSet& Foam::dimAcceleration = dimensions::acceleration;

const Foam::dimensionSet& Foam::dimDensity = dimensions::density;
const Foam::dimensionSet& Foam::dimForce = dimensions::force;
const Foam::dimensionSet& Foam::dimEnergy = dimensions::energy;
const Foam::dimensionSet& Foam::dimPower = dimensions::power;

const Foam::dimensionSet& Foam::dimPressure = dimensions::pressure;
const Foam::dimensionSet& Foam::dimSpecificHeatCapacity =
    dimensions::specificHeatCapacity;
const Foam::dimensionSet& Foam::dimKinematicViscosity =
    dimensions::kinematicViscosity;
const Foam::dimensionSet& Foam::dimDynamicViscosity =
    dimensions::dynamicViscosity;
const Foam::dimensionSet& Foam::dimThermalConductivity =
    dimensions::thermalConductivity;

const Foam::dimensionSet& Foam::dimVolumetricFlux = dimensions::volumetricFlux;
const Foam::dimensionSet& Foam::dimMassFlux = dimensions::massFlux;

// ************************************************************************* //
