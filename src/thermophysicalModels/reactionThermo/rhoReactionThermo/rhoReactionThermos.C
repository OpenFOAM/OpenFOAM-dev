/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "makeReactionThermo.H"

#include "rhoReactionThermo.H"
#include "heRhoThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "incompressiblePerfectGas.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "rhoConst.H"
#include "perfectFluid.H"
#include "adiabaticPerfectFluid.H"
#include "Boussinesq.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "homogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "veryInhomogeneousMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "singleStepReactingMixture.H"
#include "singleComponentMixture.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    homogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    inhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    veryInhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);


makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    homogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    inhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    veryInhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    incompressiblePerfectGas,
    specie
);

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    incompressiblePerfectGas,
    specie
);

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    incompressiblePerfectGas,
    specie
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Multi-component thermo for internal energy

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constGasEThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    gasEThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constIncompressibleGasEThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    incompressibleGasEThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    icoPoly8EThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constFluidEThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constAdiabaticFluidEThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constEThermoPhysics
);


// Reaction thermo for internal energy

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    reactingMixture,
    constGasEThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    reactingMixture,
    gasEThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    reactingMixture,
    constIncompressibleGasEThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    reactingMixture,
    incompressibleGasEThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    reactingMixture,
    icoPoly8EThermoPhysics
);


// Single-step reaction thermo for internal energy

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    singleStepReactingMixture,
    gasEThermoPhysics
);


// Single-component thermo for internal energy

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constGasEThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    gasEThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constIncompressibleGasEThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    incompressibleGasEThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    icoPoly8EThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constFluidEThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constAdiabaticFluidEThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constEThermoPhysics
);

makeReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    Boussinesq,
    specie
);

makeReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    Boussinesq,
    specie
);



// Multi-component thermo for sensible enthalpy

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constGasHThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    gasHThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constIncompressibleGasHThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    incompressibleGasHThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    icoPoly8HThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constFluidHThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constAdiabaticFluidHThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constHThermoPhysics
);


// Reaction thermo for sensible enthalpy

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    reactingMixture,
    constGasHThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    reactingMixture,
    gasHThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    reactingMixture,
    constIncompressibleGasHThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    reactingMixture,
    incompressibleGasHThermoPhysics
);

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    reactingMixture,
    icoPoly8HThermoPhysics
);


// Single-step reaction thermo for sensible enthalpy

makeReactionMixtureThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    singleStepReactingMixture,
    gasHThermoPhysics
);


// Single-component thermo for sensible enthalpy

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constGasHThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    gasHThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constIncompressibleGasHThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    incompressibleGasHThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    icoPoly8HThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constFluidHThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constAdiabaticFluidHThermoPhysics
);

makeReactionMixtureThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constHThermoPhysics
);

makeReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    Boussinesq,
    specie
);

makeReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    Boussinesq,
    specie
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
