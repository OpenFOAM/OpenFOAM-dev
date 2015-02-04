/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#define makeCombustionTypesThermo(CombModel, CombType, Thermo, Table)          \
                                                                               \
    typedef Foam::combustionModels::CombModel                                  \
        <Foam::combustionModels::CombType, Foam::Thermo>                       \
        CombModel##CombType##Thermo;                                           \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CombModel##CombType##Thermo,                                           \
        #CombModel"<"#CombType","#Thermo">",                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace combustionModels                                             \
        {                                                                      \
            typedef CombModel<CombType, Thermo> CombModel##CombType##Thermo;   \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                Table,                                                         \
                CombModel##CombType##Thermo,                                   \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


#define makeCombustionTypes(CombModel, CombType, Table)                        \
                                                                               \
    typedef Foam::combustionModels::CombModel                                  \
        <Foam::combustionModels::CombType>                                     \
        CombModel##CombType;                                                   \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CombModel##CombType,                                                   \
        #CombModel"<"#CombType">",                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace combustionModels                                             \
        {                                                                      \
            typedef CombModel<CombType> CombModel##CombType;                   \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                Table,                                                         \
                CombModel##CombType,                                           \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
