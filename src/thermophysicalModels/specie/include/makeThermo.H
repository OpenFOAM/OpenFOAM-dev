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

InClass
    Foam::fluidThermo

Description
    Macros for creating basic fluid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeThermo_H
#define makeThermo_H

#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineThermo(BaseThermo, CThermo, Mixture, ThermoPhys)                 \
                                                                               \
    typedef                                                                    \
        CThermo                                                                \
        <                                                                      \
            BaseThermo,                                                        \
            Mixture<ThermoPhys>                                                \
        >                                                                      \
        CThermo##Mixture##ThermoPhys;                                          \
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


#define addThermo(BaseThermo, CThermoMixtureThermoPhys)                        \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        BaseThermo,                                                            \
        CThermoMixtureThermoPhys,                                              \
        fvMesh                                                                 \
    )


#define makeThermo(BaseThermo, CThermo, Mixture, ThermoPhys)                   \
                                                                               \
    defineThermo(BaseThermo, CThermo, Mixture, ThermoPhys);                    \
                                                                               \
    addThermo(BaseThermo, CThermo##Mixture##ThermoPhys)


#define makeThermos(BaseThermo, CThermo, Mixture, ThermoPhys)                  \
                                                                               \
    defineThermo(BaseThermo, CThermo, Mixture, ThermoPhys);                    \
                                                                               \
    addThermo(basicThermo, CThermo##Mixture##ThermoPhys);                      \
    addThermo(fluidThermo, CThermo##Mixture##ThermoPhys);                      \
    addThermo(BaseThermo, CThermo##Mixture##ThermoPhys)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
