/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#ifndef forLiquids_H
#define forLiquids_H

#include "specie.H"

#include "adiabaticPerfectFluid.H"
#include "rhoConst.H"
#include "rPolynomial.H"

#include "eConstThermo.H"
#include "hConstThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "constTransport.H"

#include "thermo.H"

#include "forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forLiquidEqns(Mu, He, Cp, Macro, Args...)                              \
    forThermo(Mu, He, Cp, adiabaticPerfectFluid, specie, Macro, Args);         \
    forThermo(Mu, He, Cp, rhoConst, specie, Macro, Args);                      \
    forThermo(Mu, He, Cp, rPolynomial, specie, Macro, Args)

#define forLiquidEnergiesAndThermos(Mu, Macro, Args...)                        \
    forLiquidEqns(Mu, sensibleEnthalpy, hConstThermo, Macro, Args);            \
    forLiquidEqns(Mu, sensibleInternalEnergy, eConstThermo, Macro, Args);      \
    forLiquidEqns(Mu, sensibleInternalEnergy, hConstThermo, Macro, Args)

#define forLiquidTransports(Macro, Args...)                                    \
    forLiquidEnergiesAndThermos(constTransport, Macro, Args);

#define forLiquids(Macro, Args...)                                             \
    forLiquidTransports(Macro, Args)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forCoeffLiquidEqns(Mu, He, Cp, Macro, Args...)                         \
    forThermo(Mu, He, Cp, rhoConst, specie, Macro, Args);                      \
    forThermo(Mu, He, Cp, rPolynomial, specie, Macro, Args)

#define forCoeffLiquidEnergiesAndThermos(Mu, Macro, Args...)                   \
    forCoeffLiquidEqns(Mu, sensibleEnthalpy, hConstThermo, Macro, Args);       \
    forCoeffLiquidEqns(Mu, sensibleInternalEnergy, eConstThermo, Macro, Args); \
    forCoeffLiquidEqns(Mu, sensibleInternalEnergy, hConstThermo, Macro, Args)

#define forCoeffLiquidTransports(Macro, Args...)                               \
    forCoeffLiquidEnergiesAndThermos(constTransport, Macro, Args);

#define forCoeffLiquids(Macro, Args...)                                        \
    forCoeffLiquidTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
