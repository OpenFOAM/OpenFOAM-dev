/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#ifndef forGases_H
#define forGases_H

#include "specie.H"

#include "Boussinesq.H"
#include "incompressiblePerfectGas.H"
#include "perfectGas.H"
#include "PengRobinsonGas.H"

#include "eConstThermo.H"
#include "hConstThermo.H"
#include "janafThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "thermo.H"

#include "forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forGasEquations(Mu, He, Cp, Macro, Args...)                            \
    forThermo(Mu, He, Cp, Boussinesq, specie, Macro, Args);                    \
    forThermo(Mu, He, Cp, incompressiblePerfectGas, specie, Macro, Args);      \
    forThermo(Mu, He, Cp, PengRobinsonGas, specie, Macro, Args);               \
    forThermo(Mu, He, Cp, perfectGas, specie, Macro, Args)

#define forGasEnergiesAndThermos(Mu, Macro, Args...)                           \
    forGasEquations(Mu, sensibleEnthalpy, hConstThermo, Macro, Args);          \
    forGasEquations(Mu, sensibleEnthalpy, janafThermo, Macro, Args);           \
    forGasEquations(Mu, sensibleInternalEnergy, eConstThermo, Macro, Args);    \
    forGasEquations(Mu, sensibleInternalEnergy, hConstThermo, Macro, Args);    \
    forGasEquations(Mu, sensibleInternalEnergy, janafThermo, Macro, Args)

#define forGasTransports(Macro, Args...)                                       \
    forGasEnergiesAndThermos(constTransport, Macro, Args);                     \
    forGasEnergiesAndThermos(sutherlandTransport, Macro, Args)

#define forGases(Macro, Args...)                                               \
    forGasTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
