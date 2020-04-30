/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#ifndef makeCombustionModel_H
#define makeCombustionModel_H

#include "CombustionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineCombustionModel(ReactionThermo)                                  \
                                                                               \
    typedef CombustionModel<ReactionThermo> CombustionModel##ReactionThermo;   \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CombustionModel##ReactionThermo,                                       \
        (                                                                      \
            word(CombustionModel##ReactionThermo::typeName_()) + "<"           \
         + ReactionThermo::typeName + ">"                                      \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        CombustionModel##ReactionThermo,                                       \
        dictionary                                                             \
    );

#define makeCombustionModel(Model, ReactionThermo)                             \
                                                                               \
    typedef combustionModels::Model<ReactionThermo> Model##ReactionThermo;     \
    typedef CombustionModel<ReactionThermo> Model##ReactionThermo##Base;       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Model##ReactionThermo,                                                 \
        (                                                                      \
            word(Model##ReactionThermo::typeName_()) + "<"                     \
          + ReactionThermo::typeName + ">"                                     \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Model##ReactionThermo##Base,                                           \
        Model##ReactionThermo,                                                 \
        dictionary                                                             \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
