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

\*---------------------------------------------------------------------------*/

#ifndef makeReactionThermo_H
#define makeReactionThermo_H

#include "addToRunTimeSelectionTable.H"
#include "SpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineReactionMixtureThermo(BaseReactionThermo,CThermo,Mixture,ThermoPhys) \
                                                                               \
    typedef CThermo                                                            \
    <                                                                          \
        BaseReactionThermo,                                                    \
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
        (#CThermo"<" + Mixture<ThermoPhys>::typeName() + ">").c_str(),         \
        0                                                                      \
    )


#define addReactionMixtureThermo(BaseThermo,CThermoMixtureThermoPhys)          \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        BaseThermo,                                                            \
        CThermoMixtureThermoPhys,                                              \
        fvMesh                                                                 \
    )


#define makeReactionMixtureThermos(BaseThermo,BaseReactionThermo,CThermo,Mixture,ThermoPhys) \
                                                                               \
    defineReactionMixtureThermo(BaseReactionThermo,CThermo,Mixture,ThermoPhys);\
                                                                               \
    addReactionMixtureThermo(basicThermo, CThermo##Mixture##ThermoPhys);       \
    addReactionMixtureThermo(fluidThermo, CThermo##Mixture##ThermoPhys);       \
    addReactionMixtureThermo(BaseThermo, CThermo##Mixture##ThermoPhys);        \
    addReactionMixtureThermo(BaseReactionThermo, CThermo##Mixture##ThermoPhys)


#define makeReactionMixtureThermo(BaseReactionThermo,CThermo,Mixture,ThermoPhys) \
                                                                               \
    defineReactionMixtureThermo(BaseReactionThermo,CThermo,Mixture,ThermoPhys);\
                                                                               \
    addReactionMixtureThermo(BaseReactionThermo, CThermo##Mixture##ThermoPhys)


#define defineThermoPhys(Transport,Type,Thermo,EqnOfState,Specie)              \
                                                                               \
    typedef                                                                    \
        Transport                                                              \
        <                                                                      \
            species::thermo                                                    \
            <                                                                  \
                Thermo                                                         \
                <                                                              \
                    EqnOfState                                                 \
                    <                                                          \
                        Specie                                                 \
                    >                                                          \
                >,                                                             \
                Type                                                           \
            >                                                                  \
        >                                                                      \
        Transport##Type##Thermo##EqnOfState##Specie                            \


#define makeReactionThermos(BaseThermo,BaseReactionThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
                                                                               \
    defineThermoPhys(Transport,Type,Thermo,EqnOfState,Specie);                 \
                                                                               \
    makeReactionMixtureThermos                                                 \
    (                                                                          \
        BaseThermo,                                                            \
        BaseReactionThermo,                                                    \
        CThermo,                                                               \
        Mixture,                                                               \
        Transport##Type##Thermo##EqnOfState##Specie                            \
    )


#define makeReactionThermo(BaseReactionThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
                                                                               \
    defineThermoPhys(Transport,Type,Thermo,EqnOfState,Specie);                 \
                                                                               \
    makeReactionMixtureThermo                                                  \
    (                                                                          \
        BaseReactionThermo,                                                    \
        CThermo,                                                               \
        Mixture,                                                               \
        Transport##Type##Thermo##EqnOfState##Specie                            \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
