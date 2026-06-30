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

Global
    dimensionSets

Description
    Predefined dimensions

SourceFiles
    dimensions.C

\*---------------------------------------------------------------------------*/

#ifndef dimensions_H
#define dimensions_H

#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dimensions
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern const dimensionSet invalid;

extern const dimensionSet dimless;

extern const dimensionSet mass;
extern const dimensionSet length;
extern const dimensionSet time;
extern const dimensionSet temperature;
extern const dimensionSet moles;
extern const dimensionSet current;
extern const dimensionSet luminousIntensity;

extern const dimensionSet area;
extern const dimensionSet volume;

extern const dimensionSet rate;

extern const dimensionSet velocity;
extern const dimensionSet momentum;
extern const dimensionSet acceleration;

extern const dimensionSet density;
extern const dimensionSet force;
extern const dimensionSet energy;
extern const dimensionSet power;

extern const dimensionSet pressure;
extern const dimensionSet kinematicPressure;
extern const dimensionSet compressibility;

extern const dimensionSet gasConstant;
extern const dimensionSet specificHeatCapacity;

extern const dimensionSet kinematicViscosity;
extern const dimensionSet dynamicViscosity;
extern const dimensionSet thermalConductivity;
extern const dimensionSet dynamicDiffusivity;

extern const dimensionSet turbulentKineticEnergy;
extern const dimensionSet ReynoldsStress;
extern const dimensionSet turbulentEpsilon;
extern const dimensionSet turbulentOmega;
extern const dimensionSet turbulentViscosity;

extern const dimensionSet volumetricFlux;
extern const dimensionSet massFlux;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Table of dimensions
extern const HashTable<dimensionSet> table;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dimensions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using dimensions::dimless;

extern const dimensionSet& dimMass;
extern const dimensionSet& dimLength;
extern const dimensionSet& dimTime;
extern const dimensionSet& dimTemperature;
extern const dimensionSet& dimMoles;
extern const dimensionSet& dimCurrent;
extern const dimensionSet& dimLuminousIntensity;

extern const dimensionSet& dimArea;
extern const dimensionSet& dimVolume;

extern const dimensionSet& dimRate;

extern const dimensionSet& dimVelocity;
extern const dimensionSet& dimMomentum;
extern const dimensionSet& dimAcceleration;

extern const dimensionSet& dimDensity;
extern const dimensionSet& dimForce;
extern const dimensionSet& dimEnergy;
extern const dimensionSet& dimPower;

extern const dimensionSet& dimPressure;
extern const dimensionSet& dimKinematicPressure;
extern const dimensionSet& dimCompressibility;
extern const dimensionSet& dimGasConstant;
extern const dimensionSet& dimSpecificHeatCapacity;
extern const dimensionSet& dimKinematicViscosity;
extern const dimensionSet& dimDynamicViscosity;
extern const dimensionSet& dimThermalConductivity;

extern const dimensionSet& dimVolumetricFlux;
extern const dimensionSet& dimMassFlux;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
