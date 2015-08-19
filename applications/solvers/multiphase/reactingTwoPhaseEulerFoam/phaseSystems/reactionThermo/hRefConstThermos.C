/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
#include "makeThermo.H"

#include "rhoReactionThermo.H"
#include "heRhoThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "perfectFluid.H"
#include "rhoConst.H"
#include "incompressiblePerfectGas.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "absoluteEnthalpy.H"
#include "absoluteInternalEnergy.H"
#include "thermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "pureMixture.H"
#include "homogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "veryInhomogeneousMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "singleStepReactingMixture.H"

#include "thermoPhysicsTypes.H"

#include "hRefConstThermo.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Thermo type typedefs:

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectGas<specie>
        >,
        sensibleEnthalpy
    >
> constRefGasHThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectFluid<specie>
        >,
        sensibleEnthalpy
    >
> constRefFluidHThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectGas<specie>
        >,
        sensibleInternalEnergy
    >
> constRefGasEThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectFluid<specie>
        >,
        sensibleInternalEnergy
    >
> constRefFluidEThermoPhysics;


// pureMixture, sensibleEnthalpy:

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectFluid,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    rhoConst,
    specie
);


// pureMixture, sensibleInternalEnergy:

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hRefConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hRefConstThermo,
    perfectFluid,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hRefConstThermo,
    rhoConst,
    specie
);


// multiComponentMixture, sensibleInternalEnergy:

makeReactionMixtureThermo
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefGasEThermoPhysics
);

makeReactionMixtureThermo
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefFluidEThermoPhysics
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
