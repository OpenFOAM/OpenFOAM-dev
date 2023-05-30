/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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
    Macros for creating basic thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeThermo_H
#define makeThermo_H

#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineThermo(BaseThermo, Mixture, ThermoPhysics)                       \
                                                                               \
    typedef                                                                    \
        namedHeThermo<BaseThermo::heThermoType<Mixture<ThermoPhysics>>>        \
        BaseThermo##HeThermo##Mixture##ThermoPhysics;                          \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        BaseThermo##HeThermo##Mixture##ThermoPhysics,                          \
        (                                                                      \
            BaseThermo::heThermoName + "<"                                     \
          + Mixture<ThermoPhysics>::typeName() + ">"                           \
        ).c_str(),                                                             \
        0                                                                      \
    );


#define addThermo(SelectThermo, BaseThermo, Mixture, ThermoPhysics)            \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        SelectThermo,                                                          \
        BaseThermo##HeThermo##Mixture##ThermoPhysics,                          \
        fvMesh                                                                 \
    );


#define makeThermo(BaseThermo, Mixture, ThermoPhysics)                         \
                                                                               \
    defineThermo(BaseThermo, Mixture, ThermoPhysics);                          \
                                                                               \
    addThermo(basicThermo, BaseThermo, Mixture, ThermoPhysics);                \
    addThermo(BaseThermo, BaseThermo, Mixture, ThermoPhysics);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
