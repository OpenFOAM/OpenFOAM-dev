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

const Foam::dimensionSet Foam::dimensions::mass(1, 0, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimensions::length(0, 1, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimensions::time(0, 0, 1, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimensions::temperature(0, 0, 0, 1, 0, 0, 0);
const Foam::dimensionSet Foam::dimensions::moles(0, 0, 0, 0, 1, 0, 0);
const Foam::dimensionSet Foam::dimensions::current(0, 0, 0, 0, 0, 1, 0);
const Foam::dimensionSet Foam::dimensions::luminousIntensity
(
    0, 0, 0, 0, 0, 0, 1
);

const Foam::dimensionSet Foam::dimensions::area(sqr(length));
const Foam::dimensionSet Foam::dimensions::volume(pow3(length));

const Foam::dimensionSet Foam::dimensions::rate(dimless/time);

const Foam::dimensionSet Foam::dimensions::velocity(length/time);
const Foam::dimensionSet Foam::dimensions::momentum(mass*velocity);
const Foam::dimensionSet Foam::dimensions::acceleration(velocity/time);

const Foam::dimensionSet Foam::dimensions::density(mass/volume);
const Foam::dimensionSet Foam::dimensions::momentumDensity(momentum/volume);
const Foam::dimensionSet Foam::dimensions::force(mass*acceleration);
const Foam::dimensionSet Foam::dimensions::energy(force*length);
const Foam::dimensionSet Foam::dimensions::kineticEnergy(energy);
const Foam::dimensionSet Foam::dimensions::kineticEnergyDensity
(
    kineticEnergy/volume
);
const Foam::dimensionSet Foam::dimensions::power(energy/time);

const Foam::dimensionSet Foam::dimensions::pressure(force/area);
const Foam::dimensionSet Foam::dimensions::kinematicPressure(pressure/density);
const Foam::dimensionSet Foam::dimensions::compressibility(density/pressure);

const Foam::dimensionSet Foam::dimensions::gasConstant(energy/mass/temperature);
const Foam::dimensionSet Foam::dimensions::specificHeatCapacity(gasConstant);

const Foam::dimensionSet Foam::dimensions::kinematicViscosity(area/time);
const Foam::dimensionSet Foam::dimensions::dynamicViscosity
(
    density*kinematicViscosity
);
const Foam::dimensionSet Foam::dimensions::thermalConductivity
(
    power/length/temperature
);
const Foam::dimensionSet Foam::dimensions::dynamicDiffusivity
(
    dynamicViscosity
);

const Foam::dimensionSet Foam::dimensions::turbulentKineticEnergy
(
    sqr(velocity)
);
const Foam::dimensionSet Foam::dimensions::ReynoldsStress
(
    turbulentKineticEnergy
);
const Foam::dimensionSet Foam::dimensions::turbulentEpsilon
(
    turbulentKineticEnergy/time
);
const Foam::dimensionSet Foam::dimensions::turbulentOmega
(
    rate
);
const Foam::dimensionSet Foam::dimensions::turbulentViscosity
(
    kinematicViscosity
);

const Foam::dimensionSet Foam::dimensions::volumetricFlux(area*velocity);
const Foam::dimensionSet Foam::dimensions::massFlux(density*volumetricFlux);
const Foam::dimensionSet Foam::dimensions::heatFlux(power/area);

const Foam::dimensionSet Foam::dimensions::magneticFluxDensity
(
    force/(length*current)
);
const Foam::dimensionSet Foam::dimensions::magneticFluxPressure
(
    magneticFluxDensity*velocity
);

const Foam::HashTable<Foam::dimensionSet> Foam::dimensions::table
{
    {
        Foam::dimensionSet::dimensionTypeNames_[Foam::dimensionSet::MASS],
        Foam::dimensions::mass
    },
    {
        Foam::dimensionSet::dimensionTypeNames_[Foam::dimensionSet::LENGTH],
        Foam::dimensions::length
    },
    {
        Foam::dimensionSet::dimensionTypeNames_[Foam::dimensionSet::TIME],
        Foam::dimensions::time
    },
    {
        Foam::dimensionSet::dimensionTypeNames_
        [
            Foam::dimensionSet::TEMPERATURE
        ],
        Foam::dimensions::temperature
    },
    {
        Foam::dimensionSet::dimensionTypeNames_[Foam::dimensionSet::MOLES],
        Foam::dimensions::moles
    },
    {
        Foam::dimensionSet::dimensionTypeNames_[Foam::dimensionSet::CURRENT],
        Foam::dimensions::current
    },
    {
        Foam::dimensionSet::dimensionTypeNames_
        [
            Foam::dimensionSet::LUMINOUS_INTENSITY
        ],
        Foam::dimensions::luminousIntensity
    },
    {"area", Foam::dimensions::area},
    {"volume", Foam::dimensions::volume},

    {"rate", Foam::dimensions::rate},

    {"velocity", Foam::dimensions::velocity},
    {"momentum", Foam::dimensions::momentum},
    {"acceleration", Foam::dimensions::acceleration},

    {"density", Foam::dimensions::density},
    {"momentumDensity", Foam::dimensions::momentumDensity},
    {"force", Foam::dimensions::force},
    {"energy", Foam::dimensions::energy},
    {"kineticEnergy", Foam::dimensions::kineticEnergy},
    {"kineticEnergyDensity", Foam::dimensions::kineticEnergyDensity},
    {"power", Foam::dimensions::power},

    {"pressure", Foam::dimensions::pressure},
    {"kinematicPressure", Foam::dimensions::kinematicPressure},
    {"compressibility", Foam::dimensions::compressibility},

    {"gasConstant", Foam::dimensions::gasConstant},
    {"specificHeatCapacity", Foam::dimensions::specificHeatCapacity},

    {"kinematicViscosity", Foam::dimensions::kinematicViscosity},
    {"dynamicViscosity", Foam::dimensions::dynamicViscosity},
    {"thermalConductivity", Foam::dimensions::thermalConductivity},
    {"dynamicDiffusivity", Foam::dimensions::dynamicDiffusivity},

    {"turbulentKineticEnergy", Foam::dimensions::turbulentKineticEnergy},
    {"ReynoldsStress", Foam::dimensions::ReynoldsStress},
    {"turbulentEpsilon", Foam::dimensions::turbulentEpsilon},
    {"turbulentOmega", Foam::dimensions::turbulentOmega},
    {"turbulentViscosity", Foam::dimensions::turbulentViscosity},

    {"volumetricFlux", Foam::dimensions::volumetricFlux},
    {"massFlux", Foam::dimensions::massFlux},
    {"heatFlux", Foam::dimensions::heatFlux},

    {"magneticFluxDensity", Foam::dimensions::magneticFluxDensity},
    {"magneticFluxPressure", Foam::dimensions::magneticFluxPressure}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::dimensionSet& Foam::dimMass = dimensions::mass;
const Foam::dimensionSet& Foam::dimLength = dimensions::length;
const Foam::dimensionSet& Foam::dimTime = dimensions::time;
const Foam::dimensionSet& Foam::dimTemperature = dimensions::temperature;
const Foam::dimensionSet& Foam::dimMoles = dimensions::moles;
const Foam::dimensionSet& Foam::dimCurrent = dimensions::current;
const Foam::dimensionSet& Foam::dimLuminousIntensity =
    dimensions::luminousIntensity;

const Foam::dimensionSet& Foam::dimArea = dimensions::area;
const Foam::dimensionSet& Foam::dimVolume = dimensions::volume;

const Foam::dimensionSet& Foam::dimRate = dimensions::rate;

const Foam::dimensionSet& Foam::dimVelocity = dimensions::velocity;
const Foam::dimensionSet& Foam::dimMomentum = dimensions::momentum;
const Foam::dimensionSet& Foam::dimAcceleration = dimensions::acceleration;

const Foam::dimensionSet& Foam::dimDensity = dimensions::density;
const Foam::dimensionSet& Foam::dimForce = dimensions::force;
const Foam::dimensionSet& Foam::dimEnergy = dimensions::energy;
const Foam::dimensionSet& Foam::dimPower = dimensions::power;

const Foam::dimensionSet& Foam::dimPressure = dimensions::pressure;
const Foam::dimensionSet& Foam::dimKinematicPressure =
    dimensions::kinematicPressure;
const Foam::dimensionSet& Foam::dimCompressibility =
    dimensions::compressibility;
const Foam::dimensionSet& Foam::dimGasConstant = dimensions::gasConstant;
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
