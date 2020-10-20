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

#define defineChemistryTabulationMethod(nullArg, ThermoPhysics)                \
                                                                               \
    typedef chemistryTabulationMethod<ThermoPhysics>                           \
        chemistryTabulationMethod##ThermoPhysics;                              \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryTabulationMethod##ThermoPhysics,                              \
        (                                                                      \
            word(chemistryTabulationMethod##ThermoPhysics::typeName_())        \
          + "<" + ThermoPhysics::typeName() + ">"                              \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        chemistryTabulationMethod##ThermoPhysics,                              \
        dictionary                                                             \
    )

#define makeChemistryTabulationMethod(Method, ThermoPhysics)                   \
                                                                               \
    typedef chemistryTabulationMethod<ThermoPhysics>                           \
        chemistryTabulationMethod##ThermoPhysics;                              \
                                                                               \
    typedef chemistryTabulationMethods::Method<ThermoPhysics>                  \
        chemistryTabulationMethod##Method##ThermoPhysics;                      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryTabulationMethod##Method##ThermoPhysics,                      \
        (                                                                      \
            word(chemistryTabulationMethod##Method##ThermoPhysics::typeName_())\
         + "<" + ThermoPhysics::typeName() + ">"                               \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        chemistryTabulationMethod##ThermoPhysics,                              \
        chemistryTabulationMethod##Method##ThermoPhysics,                      \
        dictionary                                                             \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
