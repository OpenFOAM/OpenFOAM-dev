/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#ifndef makeReactionThermo_H
#define makeReactionThermo_H

#include "makeThermo.H"
#include "SpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineReactionThermo(BaseReactionThermo, CThermo, Mixture, ThermoPhys) \
                                                                               \
    typedef CThermo                                                            \
    <                                                                          \
        BaseReactionThermo::composite,                                         \
        SpecieMixture                                                          \
        <                                                                      \
            Mixture                                                            \
            <                                                                  \
                ThermoPhys                                                     \
            >                                                                  \
        >                                                                      \
    > CThermo##Mixture##ThermoPhys;                                            \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CThermo##Mixture##ThermoPhys,                                          \
        (                                                                      \
            word(CThermo##Mixture##ThermoPhys::typeName_()) + "<"              \
          + Mixture<ThermoPhys>::typeName() + ">"                              \
        ).c_str(),                                                             \
        0                                                                      \
    )


#define addReactionThermo(BaseReactionThermo, CThermo, Mixture, ThermoPhys)    \
                                                                               \
    addThermo(BaseReactionThermo, CThermo, Mixture, ThermoPhys)


#define makeReactionThermo(BaseReactionThermo, CThermo, Mixture, ThermoPhys)   \
                                                                               \
    defineReactionThermo(BaseReactionThermo, CThermo, Mixture, ThermoPhys);    \
                                                                               \
    addReactionThermo(fluidReactionThermo, CThermo, Mixture, ThermoPhys);      \
    addReactionThermo(BaseReactionThermo, CThermo, Mixture, ThermoPhys)


#define makeReactionThermos(                                                   \
    BaseThermo, BaseReactionThermo, CThermo, Mixture, ThermoPhys)              \
                                                                               \
    defineReactionThermo(BaseReactionThermo, CThermo, Mixture, ThermoPhys);    \
                                                                               \
    addReactionThermo(basicThermo, CThermo, Mixture, ThermoPhys);              \
    addReactionThermo(fluidThermo, CThermo, Mixture, ThermoPhys);              \
    addReactionThermo(BaseThermo, CThermo, Mixture, ThermoPhys);               \
    addReactionThermo(fluidReactionThermo, CThermo, Mixture, ThermoPhys);      \
    addReactionThermo(BaseReactionThermo, CThermo, Mixture, ThermoPhys)


#define makePsiuReactionThermo(Mixture, ThermoPhys)                            \
                                                                               \
    defineReactionThermo                                                       \
    (                                                                          \
        psiuReactionThermo,                                                    \
        heheuPsiThermo,                                                        \
        Mixture,                                                               \
        ThermoPhys                                                             \
    );                                                                         \
                                                                               \
    addThermo(psiuReactionThermo, heheuPsiThermo, Mixture, ThermoPhys)


#define makePsiuReactionThermos(Mixture, ThermoPhys)                           \
                                                                               \
    defineReactionThermo                                                       \
    (                                                                          \
        psiuReactionThermo,                                                    \
        heheuPsiThermo,                                                        \
        Mixture,                                                               \
        ThermoPhys                                                             \
    );                                                                         \
                                                                               \
    addThermo(basicThermo, heheuPsiThermo, Mixture, ThermoPhys);               \
    addThermo(fluidThermo, heheuPsiThermo, Mixture, ThermoPhys);               \
    addThermo(psiThermo, heheuPsiThermo, Mixture, ThermoPhys);                 \
    addThermo(psiuReactionThermo, heheuPsiThermo, Mixture, ThermoPhys)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
