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
    Foam::makeSolidReaction

Description
    Macros for instantiating reactions for solid

\*---------------------------------------------------------------------------*/

#ifndef makeSolidReaction_H
#define makeSolidReaction_H

#include "solidReaction.H"
#include "IrreversibleReaction.H"
#include "Reaction.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSolidReaction(ReactionType, Thermo, ReactionRate)                  \
                                                                               \
    typedef solidReaction<Thermo> solidReaction##Thermo;                       \
                                                                               \
    typedef Reaction<Thermo> Reaction##Thermo;                                 \
                                                                               \
    typedef ReactionType<solidReaction, Thermo, ReactionRate>                  \
        ReactionType##Thermo##ReactionRate;                                    \
                                                                               \
    defineTemplateRunTimeSelectionTable(Reaction##Thermo, dictionary);         \
                                                                               \
    defineTemplateTypeNameAndDebug(solidReaction##Thermo, 0);                  \
    defineTemplateTypeNameAndDebug(Reaction##Thermo, 0);                       \
                                                                               \
    template<>                                                                 \
    const word ReactionType##Thermo##ReactionRate::typeName                    \
    (                                                                          \
        ReactionType::typeName_()                                              \
      + ReactionRate::type()                                                   \
      + solidReaction##Thermo::typeName_()                                     \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Reaction##Thermo,                                                      \
        ReactionType##Thermo##ReactionRate,                                    \
        dictionary                                                             \
    );

#define makeSolidIRReactions(Thermo, ReactionRate)                             \
                                                                               \
    makeSolidReaction(IrreversibleReaction, Thermo, ReactionRate)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
