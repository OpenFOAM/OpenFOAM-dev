/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::makeReaction

Description
    Macros for instantiating reactions on given thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeReaction_H
#define makeReaction_H

#include "Reaction.H"

#include "IrreversibleReaction.H"
#include "ReversibleReaction.H"
#include "NonEquilibriumReversibleReaction.H"

#include "thermo.H"

#include "sutherlandTransport.H"
#include "janafThermo.H"
#include "perfectGas.H"

#include "polynomialTransport.H"
#include "hPolynomialThermo.H"
#include "icoPolynomial.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReaction(Thermo, ReactionType, ReactionRate)                       \
                                                                               \
    typedef Reaction<Thermo> Reaction##Thermo;                                 \
                                                                               \
    typedef ReactionType<Reaction, Thermo, ReactionRate>                       \
        ReactionType##Thermo##ReactionRate;                                    \
                                                                               \
    template<>                                                                 \
    const word ReactionType##Thermo##ReactionRate::typeName                    \
    (                                                                          \
        ReactionType::typeName_()                                              \
      + ReactionRate::type()                                                   \
      + Reaction##Thermo::typeName_()                                          \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Reaction##Thermo,                                                      \
        ReactionType##Thermo##ReactionRate,                                    \
        dictionary                                                             \
    );


#define makePressureDependentReaction\
(                                                                              \
    Thermo,                                                                    \
    Reaction,                                                                  \
    PressureDependentReactionRate,                                             \
    ReactionRate,                                                              \
    FallOffFunction                                                            \
)                                                                              \
                                                                               \
    typedef PressureDependentReactionRate<ReactionRate, FallOffFunction>       \
        PressureDependentReactionRate##ReactionRate##FallOffFunction;          \
                                                                               \
    makeReaction                                                               \
    (                                                                          \
        Thermo,                                                                \
        Reaction,                                                              \
        PressureDependentReactionRate##ReactionRate##FallOffFunction           \
    )


#define makeIRReactions(Thermo, ReactionRate)                                  \
                                                                               \
    makeReaction(Thermo, IrreversibleReaction, ReactionRate)                   \
                                                                               \
    makeReaction(Thermo, ReversibleReaction, ReactionRate)


#define makeIRNReactions(Thermo, ReactionRate)                                 \
                                                                               \
    makeIRReactions(Thermo, ReactionRate)                                      \
                                                                               \
    makeReaction(Thermo, NonEquilibriumReversibleReaction, ReactionRate)


#define makePressureDependentReactions(Thermo, ReactionRate, FallOffFunction)  \
                                                                               \
    makePressureDependentReaction                                              \
    (                                                                          \
        Thermo,                                                                \
        IrreversibleReaction,                                                  \
        FallOffReactionRate,                                                   \
        ReactionRate,                                                          \
        FallOffFunction                                                        \
    )                                                                          \
                                                                               \
    makePressureDependentReaction                                              \
    (                                                                          \
        Thermo,                                                                \
        ReversibleReaction,                                                    \
        FallOffReactionRate,                                                   \
        ReactionRate,                                                          \
        FallOffFunction                                                        \
    )                                                                          \
                                                                               \
    makePressureDependentReaction                                              \
    (                                                                          \
        Thermo,                                                                \
        IrreversibleReaction,                                                  \
        ChemicallyActivatedReactionRate,                                       \
        ReactionRate,                                                          \
        FallOffFunction                                                        \
    )                                                                          \
                                                                               \
    makePressureDependentReaction                                              \
    (                                                                          \
        Thermo,                                                                \
        ReversibleReaction,                                                    \
        ChemicallyActivatedReactionRate,                                       \
        ReactionRate,                                                          \
        FallOffFunction                                                        \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
