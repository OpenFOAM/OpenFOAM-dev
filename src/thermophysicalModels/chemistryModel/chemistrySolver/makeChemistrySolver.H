/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#ifndef makeChemistrySolver_H
#define makeChemistrySolver_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistrySolver(Model, ReactionThermo, ThermoPhysics)            \
                                                                               \
    typedef Model<ReactionThermo, ThermoPhysics>                               \
        Model##ReactionThermo##ThermoPhysics;                                  \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Model##ReactionThermo##ThermoPhysics,                                  \
        (                                                                      \
            word(Model##ReactionThermo##ThermoPhysics::typeName_()) + "<"      \
          + ReactionThermo::typeName + "," + ThermoPhysics::typeName() + ">"   \
        ).c_str(),                                                             \
        0                                                                      \
    );


#define makeChemistrySolver(Solver, Model, ReactionThermo, ThermoPhysics)      \
                                                                               \
    typedef Solver<Model<ReactionThermo, ThermoPhysics>>                       \
        Solver##Model##ReactionThermo##ThermoPhysics;                          \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Solver##Model##ReactionThermo##ThermoPhysics,                          \
        (                                                                      \
            word(Solver##Model##ReactionThermo##ThermoPhysics::typeName_())    \
          + "<" + word(Model<ReactionThermo, ThermoPhysics>::typeName_())      \
          + "<" + word(ReactionThermo::typeName_()) + ","                      \
          + ThermoPhysics::typeName() + ">>"                                   \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<ReactionThermo>::                                      \
        addthermoConstructorToTable                                            \
        <Solver##Model##ReactionThermo##ThermoPhysics>                         \
        add##Solver##Model##ReactionThermo##ThermoPhysics##\
thermoConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
