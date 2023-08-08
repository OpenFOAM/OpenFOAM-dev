/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Useful dimension sets

SourceFiles
    dimensionSets.C

\*---------------------------------------------------------------------------*/

#ifndef dimensionSets_H
#define dimensionSets_H

#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern const dimensionSet dimless;

extern const dimensionSet dimMass;
extern const dimensionSet dimLength;
extern const dimensionSet dimTime;
extern const dimensionSet dimTemperature;
extern const dimensionSet dimMoles;
extern const dimensionSet dimCurrent;
extern const dimensionSet dimLuminousIntensity;

extern const dimensionSet dimArea;
extern const dimensionSet dimVolume;
extern const dimensionSet dimVol;

extern const dimensionSet dimVelocity;
extern const dimensionSet dimMomentum;
extern const dimensionSet dimAcceleration;

extern const dimensionSet dimDensity;
extern const dimensionSet dimForce;
extern const dimensionSet dimEnergy;
extern const dimensionSet dimPower;

extern const dimensionSet dimPressure;
extern const dimensionSet dimCompressibility;
extern const dimensionSet dimGasConstant;
extern const dimensionSet dimSpecificHeatCapacity;
extern const dimensionSet dimViscosity;
extern const dimensionSet dimDynamicViscosity;

extern const dimensionSet dimFlux;
extern const dimensionSet dimMassFlux;


//- Add a unit
void addUnit(const dimensionedScalar& unit);


//- Set of all units
const HashTable<dimensionedScalar>& unitSet();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
