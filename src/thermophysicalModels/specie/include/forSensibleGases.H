/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#ifndef forSensibleGases_H
#define forSensibleGases_H

#include "specie.H"

#include "perfectGas.H"

#include "eConstThermo.H"
#include "hConstThermo.H"
#include "janafThermo.H"

#include "sensibleEnthalpy.H"

#include "constTransport.H"
#include "sutherlandTransport.H"
#include "polynomialTransport.H"

#include "thermo.H"

#include "forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forSensibleGasEqns(Mu, He, Cp, Macro, Args...)                         \
    forThermo(Mu, He, Cp, perfectGas, specie, Macro, Args)

#define forSensibleGasEnergiesAndThermos(Mu, Macro, Args...)                   \
    forSensibleGasEqns(Mu, sensibleEnthalpy, hConstThermo, Macro, Args);       \
    forSensibleGasEqns(Mu, sensibleEnthalpy, janafThermo, Macro, Args);

#define forSensibleGasTransports(Macro, Args...)                               \
    forSensibleGasEnergiesAndThermos(constTransport, Macro, Args);             \
    forSensibleGasEnergiesAndThermos(sutherlandTransport, Macro, Args)         \
    forSensibleGasEnergiesAndThermos(polynomialTransport, Macro, Args)

#define forSensibleGases(Macro, Args...)                                       \
    forSensibleGasTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
