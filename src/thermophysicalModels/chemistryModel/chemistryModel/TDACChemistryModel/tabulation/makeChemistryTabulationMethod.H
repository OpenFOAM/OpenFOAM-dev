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

#ifndef makeChemistryTabulationMethod_H
#define makeChemistryTabulationMethod_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistryTabulationMethod(ReactionThermo, ThermoPhysics)         \
                                                                               \
    typedef chemistryTabulationMethod<ReactionThermo, ThermoPhysics>           \
        chemistryTabulationMethod##ReactionThermo##ThermoPhysics;              \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryTabulationMethod##ReactionThermo##ThermoPhysics,              \
        (                                                                      \
            word                                                               \
            (                                                                  \
                chemistryTabulationMethod##ReactionThermo##ThermoPhysics::     \
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
        chemistryTabulationMethod##ReactionThermo##ThermoPhysics,              \
        dictionary                                                             \
    );                                                                         \

#define makeChemistryTabulationMethod(Method, ReactionThermo, ThermoPhysics)   \
                                                                               \
    typedef chemistryTabulationMethods::Method<ReactionThermo, ThermoPhysics>  \
        chemistryTabulationMethod##Method##ReactionThermo##ThermoPhysics;      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryTabulationMethod##Method##ReactionThermo##ThermoPhysics,      \
        (                                                                      \
            word                                                               \
            (                                                                  \
                chemistryTabulationMethod##Method##ReactionThermo##\
ThermoPhysics::typeName_()                                                     \
            )                                                                  \
         + "<" + word(ReactionThermo::typeName_())                             \
         + "," + ThermoPhysics::typeName() + ">"                               \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    chemistryTabulationMethod<ReactionThermo, ThermoPhysics>::                 \
        adddictionaryConstructorToTable                                        \
        <chemistryTabulationMethod##Method##ReactionThermo##ThermoPhysics>     \
        add##chemistryTabulationMethod##Method##ReactionThermo##ThermoPhysics##\
dictionaryConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
