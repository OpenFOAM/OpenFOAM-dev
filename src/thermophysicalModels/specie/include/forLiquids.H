/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "specie.H"

#include "adiabaticPerfectFluid.H"
#include "rhoConst.H"
#include "rPolynomial.H"
#include "NSRDS.H"

#include "eConstThermo.H"
#include "hConstThermo.H"
#include "NSRDSThermo.H"
#include "NSRDS14Thermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "constTransport.H"
#include "NSRDSTransport.H"

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
    forLiquidEnergiesAndThermos(constTransport, Macro, Args)

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
    forCoeffLiquidEnergiesAndThermos(constTransport, Macro, Args)

#define forCoeffLiquids(Macro, Args...)                                        \
    forCoeffLiquidTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forEnergyLiquidEqns(Mu, He, Cp, Macro, Args...)                        \
    forThermo(Mu, He, Cp, adiabaticPerfectFluid, specie, Macro, Args);         \
    forThermo(Mu, He, Cp, rhoConst, specie, Macro, Args);                      \
    forThermo(Mu, He, Cp, rPolynomial, specie, Macro, Args)

#define forEnergyLiquidEnergiesAndThermos(Mu, Macro, Args...)                  \
    forEnergyLiquidEqns(Mu, sensibleInternalEnergy, eConstThermo, Macro, Args);\
    forEnergyLiquidEqns(Mu, sensibleInternalEnergy, hConstThermo, Macro, Args)

#define forEnergyLiquidTransports(Macro, Args...)                              \
    forEnergyLiquidEnergiesAndThermos(constTransport, Macro, Args)

#define forEnergyLiquids(Macro, Args...)                                       \
    forEnergyLiquidTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forCoeffEnergyLiquidEqns(Mu, He, Cp, Macro, Args...)                   \
    forThermo(Mu, He, Cp, rhoConst, specie, Macro, Args);                      \
    forThermo(Mu, He, Cp, rPolynomial, specie, Macro, Args)

#define forCoeffEnergyLiquidEnergiesAndThermos(Mu, Args...)                    \
    forCoeffEnergyLiquidEqns(Mu, sensibleInternalEnergy, eConstThermo, Args);  \
    forCoeffEnergyLiquidEqns(Mu, sensibleInternalEnergy, hConstThermo, Args)

#define forCoeffEnergyLiquidTransports(Macro, Args...)                         \
    forCoeffEnergyLiquidEnergiesAndThermos(constTransport, Macro, Args)

#define forCoeffEnergyLiquids(Macro, Args...)                                  \
    forCoeffEnergyLiquidTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forNSRDSLiquidEqns(Mu, He, Cp, Macro, Args...)                         \
    forThermo(Mu, He, Cp, NSRDS, specie, Macro, Args);                         \

#define forNSRDSLiquidThermos(Mu, He, Macro, Args...)                          \
    forNSRDSLiquidEqns(Mu, He, NSRDSThermo, Macro, Args);                      \
    forNSRDSLiquidEqns(Mu, He, NSRDS14Thermo, Macro, Args)

#define forNSRDSLiquidEnergies(Mu, Macro, Args...)                             \
    forNSRDSLiquidThermos(Mu, sensibleEnthalpy, Macro, Args);                  \
    forNSRDSLiquidThermos(Mu, sensibleInternalEnergy, Macro, Args)

#define forNSRDSLiquidTransports(Macro, Args...)                               \
    forNSRDSLiquidEnergies(NSRDSTransport, Macro, Args)

#define forNSRDSLiquids(Macro, Args...)                                        \
    forNSRDSLiquidTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forNSRDSEnergyLiquidEqns(Mu, He, Cp, Macro, Args...)                   \
    forThermo(Mu, He, Cp, NSRDS, specie, Macro, Args);                         \

#define forNSRDSEnergyLiquidThermos(Mu, He, Macro, Args...)                    \
    forNSRDSEnergyLiquidEqns(Mu, He, NSRDSThermo, Macro, Args);                \
    forNSRDSEnergyLiquidEqns(Mu, He, NSRDS14Thermo, Macro, Args)

#define forNSRDSEnergyLiquidEnergies(Mu, Macro, Args...)                       \
    forNSRDSEnergyLiquidThermos(Mu, sensibleEnthalpy, Macro, Args);            \
    forNSRDSEnergyLiquidThermos(Mu, sensibleInternalEnergy, Macro, Args)

#define forNSRDSEnergyLiquidTransports(Macro, Args...)                         \
    forNSRDSEnergyLiquidEnergies(NSRDSTransport, Macro, Args)

#define forNSRDSEnergyLiquids(Macro, Args...)                                  \
    forNSRDSEnergyLiquidTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
