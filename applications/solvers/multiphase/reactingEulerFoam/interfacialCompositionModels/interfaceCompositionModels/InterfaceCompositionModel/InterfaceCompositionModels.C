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

namespace Foam
{
    using namespace interfaceCompositionModels;

    // Gas-side models

    // multi-component gas in the presence of a pure liquid
    makeInterfaceCompositionType
    (
        Saturated,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        gasEThermoPhysics,
        heRhoThermo,
        rhoThermo,
        pureMixture,
        constFluidEThermoPhysics
    );

    // reacting gas in the presence of a pure liquid
    makeInterfaceCompositionType
    (
        Saturated,
        heRhoThermo,
        rhoReactionThermo,
        reactingMixture,
        gasEThermoPhysics,
        heRhoThermo,
        rhoThermo,
        pureMixture,
        constFluidEThermoPhysics
    );

    // multi-component gas in the presence of a multi-component liquid
    makeSpecieInterfaceCompositionType
    (
        Saturated,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constGasEThermoPhysics,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constFluidEThermoPhysics
    );
    makeSpecieInterfaceCompositionType
    (
        NonRandomTwoLiquid,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constGasEThermoPhysics,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constFluidEThermoPhysics
    );

    // reacting gas in the presence of a multi-component liquid
    makeSpecieInterfaceCompositionType
    (
        Saturated,
        heRhoThermo,
        rhoReactionThermo,
        reactingMixture,
        constGasEThermoPhysics,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constFluidEThermoPhysics
    );
    makeSpecieInterfaceCompositionType
    (
        NonRandomTwoLiquid,
        heRhoThermo,
        rhoReactionThermo,
        reactingMixture,
        constGasEThermoPhysics,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constFluidEThermoPhysics
    );

    // Liquid-side models

    // multi-component liquid in the presence of a multi-component gas
    makeSpecieInterfaceCompositionType
    (
        Henry,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constFluidEThermoPhysics,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constGasEThermoPhysics
    );
    makeSpecieInterfaceCompositionType
    (
        Raoult,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constFluidEThermoPhysics,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constGasEThermoPhysics
    );

    // multi-component liquid in the presence of a reacting gas
    makeSpecieInterfaceCompositionType
    (
        Henry,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constFluidEThermoPhysics,
        heRhoThermo,
        rhoReactionThermo,
        reactingMixture,
        constGasEThermoPhysics
    );
    makeSpecieInterfaceCompositionType
    (
        Raoult,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constFluidEThermoPhysics,
        heRhoThermo,
        rhoReactionThermo,
        reactingMixture,
        constGasEThermoPhysics
    );
}

// ************************************************************************* //
