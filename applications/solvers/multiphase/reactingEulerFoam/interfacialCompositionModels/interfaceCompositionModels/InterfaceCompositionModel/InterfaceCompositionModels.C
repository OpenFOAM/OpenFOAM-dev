/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interfaceCompositionModel.H"
#include "InterfaceCompositionModel.H"
#include "Henry.H"
#include "NonRandomTwoLiquid.H"
#include "Raoult.H"
#include "Saturated.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "makeReactionThermo.H"

#include "thermoPhysicsTypes.H"

#include "rhoConst.H"
#include "perfectFluid.H"

#include "pureMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "SpecieMixture.H"

#include "rhoThermo.H"
#include "rhoReactionThermo.H"
#include "heRhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSpecieInterfaceCompositionModel(Model, Thermo1, Thermo2)           \
                                                                               \
    /* Composition at an interface with a multi-component mixture */           \
    makeSpecieInterfaceCompositionType                                         \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        multiComponentMixture, Thermo1,                                        \
        heRhoThermo, rhoReactionThermo,                                        \
        multiComponentMixture, Thermo2                                         \
    );                                                                         \
    makeSpecieInterfaceCompositionType                                         \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        reactingMixture, Thermo1,                                              \
        heRhoThermo, rhoReactionThermo,                                        \
        multiComponentMixture, Thermo2                                         \
    );                                                                         \
                                                                               \
    /* Composition at an interface with a reacting mixture */                  \
    makeSpecieInterfaceCompositionType                                         \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        multiComponentMixture, Thermo1,                                        \
        heRhoThermo, rhoReactionThermo,                                        \
        reactingMixture, Thermo2                                               \
    );                                                                         \
    makeSpecieInterfaceCompositionType                                         \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        reactingMixture, Thermo1,                                              \
        heRhoThermo, rhoReactionThermo,                                        \
        reactingMixture, Thermo2                                               \
    );

#define makeInterfaceCompositionModel(Model, Thermo1, Thermo2)                 \
                                                                               \
    /* Composition at an interface with a pure mixture */                      \
    makeInterfaceCompositionType                                               \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        multiComponentMixture, Thermo1,                                        \
        heRhoThermo, rhoThermo,                                                \
        pureMixture, Thermo2                                                   \
    );                                                                         \
    makeInterfaceCompositionType                                               \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        reactingMixture, Thermo1,                                              \
        heRhoThermo, rhoThermo,                                                \
        pureMixture, Thermo2                                                   \
    );                                                                         \
                                                                               \
    /* Composition at an interface with non-pure mixtures */                   \
    makeSpecieInterfaceCompositionModel(Model, Thermo1, Thermo2)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    using namespace interfaceCompositionModels;

    // Gas-side models
    makeInterfaceCompositionModel
    (
        Saturated,
        gasEThermoPhysics,
        constFluidEThermoPhysics
    );
    makeInterfaceCompositionModel
    (
        Saturated,
        constGasEThermoPhysics,
        constFluidEThermoPhysics
    );
    makeSpecieInterfaceCompositionModel
    (
        NonRandomTwoLiquid,
        gasEThermoPhysics,
        constFluidEThermoPhysics
    );
    makeSpecieInterfaceCompositionModel
    (
        NonRandomTwoLiquid,
        constGasEThermoPhysics,
        constFluidEThermoPhysics
    );

    // Liquid-side models
    makeSpecieInterfaceCompositionModel
    (
        Henry,
        constFluidEThermoPhysics,
        gasEThermoPhysics
    );
    makeSpecieInterfaceCompositionModel
    (
        Henry,
        constFluidEThermoPhysics,
        constGasEThermoPhysics
    );
    makeSpecieInterfaceCompositionModel
    (
        Raoult,
        constFluidEThermoPhysics,
        gasEThermoPhysics
    );
    makeSpecieInterfaceCompositionModel
    (
        Raoult,
        constFluidEThermoPhysics,
        constGasEThermoPhysics
    );
}

// ************************************************************************* //
