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

#ifndef makeChemistryReductionMethod_H
#define makeChemistryReductionMethod_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistryReductionMethod(ReactionThermo, ThermoPhysics)          \
                                                                               \
    typedef chemistryReductionMethod<ReactionThermo, ThermoPhysics>            \
        chemistryReductionMethod##ReactionThermo##ThermoPhysics;               \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryReductionMethod##ReactionThermo##ThermoPhysics,               \
        (                                                                      \
            word                                                               \
            (                                                                  \
                chemistryReductionMethod##ReactionThermo##ThermoPhysics::      \
                typeName_()                                                    \
            )                                                                  \
          + "<" + word(ReactionThermo::typeName_())                            \
          + "," + ThermoPhysics::typeName() + ">"                              \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        chemistryReductionMethod##ReactionThermo##ThermoPhysics,               \
        dictionary                                                             \
    );                                                                         \

#define makeChemistryReductionMethod(Method, ReactionThermo, ThermoPhysics)    \
                                                                               \
    typedef chemistryReductionMethods::Method<ReactionThermo, ThermoPhysics>   \
        chemistryReductionMethod##Method##ReactionThermo##ThermoPhysics;       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryReductionMethod##Method##ReactionThermo##ThermoPhysics,       \
        (                                                                      \
            word                                                               \
            (                                                                  \
                chemistryReductionMethod##Method##ReactionThermo##\
ThermoPhysics::typeName_()                                                     \
            )                                                                  \
         + "<" + word(ReactionThermo::typeName_())                             \
         + "," + ThermoPhysics::typeName() + ">"                               \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    chemistryReductionMethod<ReactionThermo, ThermoPhysics>::                  \
        adddictionaryConstructorToTable                                        \
        <chemistryReductionMethod##Method##ReactionThermo##ThermoPhysics>      \
        add##chemistryReductionMethod##Method##ReactionThermo##ThermoPhysics##\
dictionaryConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
