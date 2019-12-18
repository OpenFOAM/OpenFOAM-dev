/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "mixtureFraction.H"
#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeMixtureSootTypesThermo(SootModelType, ReactionThermo, Thermo)      \
                                                                               \
    typedef radiationModels::sootModels::SootModelType<ReactionThermo, Thermo> \
        SootModelType##ReactionThermo##Thermo;                                 \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SootModelType##ReactionThermo##Thermo,                                 \
        (                                                                      \
            word(SootModelType##ReactionThermo##Thermo::typeName_()) +         \
            "<"#ReactionThermo","#Thermo">"                                    \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace radiationModels                                                  \
    {                                                                          \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            sootModel,                                                         \
            SootModelType##ReactionThermo##Thermo,                             \
            dictionary                                                         \
        );                                                                     \
    }

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeMixtureSootTypesThermo
    (
        mixtureFraction,
        psiReactionThermo,
        gasHThermoPhysics
    );

    makeMixtureSootTypesThermo
    (
        mixtureFraction,
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeMixtureSootTypesThermo
    (
        mixtureFraction,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeMixtureSootTypesThermo
    (
        mixtureFraction,
        rhoReactionThermo,
        gasEThermoPhysics
    );
}

// ************************************************************************* //
