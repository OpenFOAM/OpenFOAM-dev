/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
    Foam::solidThermo

Description
    Macros for creating solid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeSolidThermo_H
#define makeSolidThermo_H

#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineSolidThermo(BaseThermo, CThermo, Mixture, ThermoPhys)            \
                                                                               \
    typedef                                                                    \
        heThermo                                                               \
        <                                                                      \
            BaseThermo::composite,                                             \
            Mixture<ThermoPhys>                                                \
        > heThermo##Mixture##ThermoPhys;                                       \
                                                                               \
    typedef                                                                    \
        CThermo                                                                \
        <                                                                      \
            BaseThermo::composite,                                             \
            Mixture<ThermoPhys>                                                \
        > CThermo##Mixture##ThermoPhys;                                        \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CThermo##Mixture##ThermoPhys,                                          \
        (                                                                      \
            #CThermo"<"#Mixture"<"                                             \
          + ThermoPhys::typeName()                                             \
          + ">>"                                                               \
        ).c_str(),                                                             \
        0                                                                      \
    )


#define addSolidThermo(BaseThermo, CThermoMixtureThermoPhys, Table)            \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        BaseThermo,                                                            \
        CThermoMixtureThermoPhys,                                              \
        Table                                                                  \
    )


#define makeSolidThermo(BaseThermo, CThermo, Mixture, ThermoPhys)              \
                                                                               \
    defineSolidThermo(BaseThermo, CThermo, Mixture, ThermoPhys);               \
                                                                               \
    addSolidThermo(BaseThermo, CThermo##Mixture##ThermoPhys, fvMesh


#define makeSolidThermos(BaseThermo, CThermo, Mixture, ThermoPhys)             \
                                                                               \
    defineSolidThermo(BaseThermo, CThermo, Mixture, ThermoPhys);               \
                                                                               \
    addSolidThermo(basicThermo, CThermo##Mixture##ThermoPhys, fvMesh);         \
    addSolidThermo(BaseThermo, CThermo##Mixture##ThermoPhys, fvMesh);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
