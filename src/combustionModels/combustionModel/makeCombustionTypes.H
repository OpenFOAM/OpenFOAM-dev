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

#ifndef makeCombustionTypes_H
#define makeCombustionTypes_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCombustion(Comp)                                                   \
                                                                               \
    typedef typename Foam::combustionModels::CombustionModel<Foam::Comp>       \
        CombustionModel##Comp;                                                 \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CombustionModel##Comp,                                                 \
        "CombustionModel<"#Comp">",                                            \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        CombustionModel##Comp,                                                 \
        dictionary                                                             \
    );


#define makeCombustionTypesThermo(CombModel, Comp, Thermo)                     \
                                                                               \
    typedef Foam::combustionModels::CombModel<Foam::Comp, Foam::Thermo>        \
        CombModel##Comp##Thermo;                                               \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CombModel##Comp##Thermo,                                               \
        (                                                                      \
            Foam::word(CombModel##Comp##Thermo::typeName_()) + "<"#Comp","     \
          + Foam::Thermo::typeName() + ">"                                     \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    Foam::combustionModels::CombustionModel<Foam::Comp>::                      \
        add##dictionary##ConstructorToTable<CombModel##Comp##Thermo>           \
        add##CombModel##Comp##Thermo##dictionary##ConstructorTo##\
CombustionModel##Comp##Table_;


#define makeCombustionTypes(CombModel, Comp)                                   \
                                                                               \
    typedef Foam::combustionModels::CombModel<Foam::Comp>                      \
        CombModel##Comp;                                                       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CombModel##Comp,                                                       \
        (Foam::word(CombModel##Comp::typeName_()) + "<"#Comp">").c_str(),      \
        0                                                                      \
    );                                                                         \
                                                                               \
    Foam::combustionModels::CombustionModel<Foam::Comp>::                      \
        add##dictionary##ConstructorToTable<CombModel##Comp>                   \
        add##CombModel##Comp##dictionary##ConstructorTo##CombustionModel##Comp\
##Table_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
