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

InClass
    Foam::makeRateReaction

Description
    Macros for instantiating reactions on given thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeRateReaction_H
#define makeRateReaction_H

#include "Irreversible_reaction.H"
#include "Reversible_reaction.H"
#include "NonEquilibriumReversible_reaction.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeRateReaction(Table, Type, Rate)                                    \
                                                                               \
    typedef reactions::Type<reactionRates::Rate> Type##Rate;                   \
                                                                               \
    template<>                                                                 \
    const word Type##Rate::typeName =                                          \
        word(Type##Rate::typeName_()) + '<'                                    \
      + reactionRates::Rate::type() + '>';                                     \
                                                                               \
    addToRunTimeSelectionTable(reaction, Type##Rate, Table)


#define makeThermoRateReaction(Table, Type, Rate, ThermoPhysics)               \
                                                                               \
    typedef reactions::Type<reactionRates::Rate, ThermoPhysics>                \
        Type##Rate##ThermoPhysics;                                             \
                                                                               \
    template<>                                                                 \
    const word Type##Rate##ThermoPhysics::typeName =                           \
        word(Type##Rate##ThermoPhysics::typeName_()) + '<'                     \
      + reactionRates::Rate::type() + ',' + ThermoPhysics::typeName() + '>';   \
                                                                               \
    addToRunTimeSelectionTable(reaction, Type##Rate##ThermoPhysics, Table)


#define makeIrreversibleReaction(Table, Rate)                                  \
    makeRateReaction(Table, Irreversible, Rate)


#define makeReversibleReaction(Table, Rate, ThermoPhysics)                     \
    makeThermoRateReaction(Table, Reversible, Rate, ThermoPhysics)


#define makeNonEquilibriumReversibleReaction(Table, Rate)                      \
    makeRateReaction(Table, NonEquilibriumReversible, Rate)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
