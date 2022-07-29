/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#ifndef makeMulticomponentThermo_H
#define makeMulticomponentThermo_H

#include "makeThermo.H"
#include "SpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineMulticomponentThermo(                                            \
    BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)                    \
                                                                               \
    typedef CThermo                                                            \
    <                                                                          \
        BaseMulticomponentThermo::composite,                                   \
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


#define addMulticomponentThermo(                                               \
    BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)                    \
                                                                               \
    addThermo(BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)


#define makeMulticomponentThermo(                                              \
    BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)                    \
                                                                               \
    defineMulticomponentThermo                                                 \
    (BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys);                  \
                                                                               \
    addMulticomponentThermo                                                    \
    (fluidMulticomponentThermo, CThermo, Mixture, ThermoPhys);                 \
                                                                               \
    addMulticomponentThermo                                                    \
    (BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)


#define makeMulticomponentThermos(                                             \
    BaseThermo, BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)        \
                                                                               \
    defineMulticomponentThermo                                                 \
    (BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys);                  \
                                                                               \
    addMulticomponentThermo(basicThermo, CThermo, Mixture, ThermoPhys);        \
    addMulticomponentThermo(fluidThermo, CThermo, Mixture, ThermoPhys);        \
    addMulticomponentThermo(BaseThermo, CThermo, Mixture, ThermoPhys);         \
                                                                               \
    addMulticomponentThermo                                                    \
    (fluidMulticomponentThermo, CThermo, Mixture, ThermoPhys);                 \
                                                                               \
    addMulticomponentThermo                                                    \
    (BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)


#define makePsiuMulticomponentThermo(Mixture, ThermoPhys)                      \
                                                                               \
    defineMulticomponentThermo                                                 \
    (                                                                          \
        psiuMulticomponentThermo,                                              \
        heheuPsiThermo,                                                        \
        Mixture,                                                               \
        ThermoPhys                                                             \
    );                                                                         \
                                                                               \
    addThermo(psiuMulticomponentThermo, heheuPsiThermo, Mixture, ThermoPhys)


#define makePsiuMulticomponentThermos(Mixture, ThermoPhys)                     \
                                                                               \
    defineMulticomponentThermo                                                 \
    (                                                                          \
        psiuMulticomponentThermo,                                              \
        heheuPsiThermo,                                                        \
        Mixture,                                                               \
        ThermoPhys                                                             \
    );                                                                         \
                                                                               \
    addThermo(basicThermo, heheuPsiThermo, Mixture, ThermoPhys);               \
    addThermo(fluidThermo, heheuPsiThermo, Mixture, ThermoPhys);               \
    addThermo(psiThermo, heheuPsiThermo, Mixture, ThermoPhys);                 \
    addThermo(psiuMulticomponentThermo, heheuPsiThermo, Mixture, ThermoPhys)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
