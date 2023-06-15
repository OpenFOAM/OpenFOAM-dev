/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

Description
    Macros for easy insertion into run-time selection tables

\*---------------------------------------------------------------------------*/

#ifndef addToRunTimeSelectionTable_H
#define addToRunTimeSelectionTable_H

#include "macros.H"

namespace Foam
{
    class word;

    bool addToRunTimeSelectionTableToC
    (
        const char* baseType,
        const char* baseTypeName,
        const word& thisTypeName,
        const char* libName
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to hash-table of functions with typename as the key
#define addToRunTimeSelectionTable\
(baseType,thisType,argNames)                                                   \
                                                                               \
    /* Add the thisType constructor function to the table */                   \
    baseType::add##argNames##ConstructorToTable<thisType>                      \
        add##thisType##argNames##ConstructorTo##baseType##Table_;              \
                                                                               \
    bool add##thisType##argNames##To##baseType##ToC_ =                         \
        addToRunTimeSelectionTableToC                                          \
        (#baseType, baseType::typeName_(), thisType::typeName, STR(LIB_NAME))


//- Add to hash-table of functions with 'lookup' as the key
#define addNamedToRunTimeSelectionTable\
(baseType,thisType,argNames,lookup)                                            \
                                                                               \
    /* Add the thisType constructor function to the table, find by lookup */   \
    baseType::add##argNames##ConstructorToTable<thisType>                      \
        add_##lookup##_##thisType##argNames##ConstructorTo##baseType##Table_   \
        (#lookup);                                                             \
                                                                               \
    bool add_##lookup##_##thisType##argNames##To##baseType##ToC_ =             \
        addToRunTimeSelectionTableToC                                          \
        (#baseType, baseType::typeName_(), #lookup, STR(LIB_NAME))


//- Add backward-compatibility duplicate of thisType, replacing oldType,
//  to hash-table of functions with 'lookup' as the key
#define addBackwardCompatibleToRunTimeSelectionTable\
(baseType,thisType,argNames,oldType,lookup)                                    \
                                                                               \
    /* Add the thisType constructor function to the table, find by lookup */   \
    baseType::add##argNames##ConstructorToTable<thisType>                      \
        add_##oldType##_##thisType##argNames##ConstructorTo##baseType##Table_  \
        (lookup);


//- Add to hash-table of functions with typename as the key
#define addRemovableToRunTimeSelectionTable\
(baseType,thisType,argNames)                                                   \
                                                                               \
    /* Add the thisType constructor function to the table */                   \
    baseType::addRemovable##argNames##ConstructorToTable<thisType>             \
        addRemovable##thisType##argNames##ConstructorTo##baseType##Table_;     \
                                                                               \
    bool addRemovable##thisType##argNames##To##baseType##ToC_ =                \
        addToRunTimeSelectionTableToC                                          \
        (#baseType, baseType::typeName_(), thisType::typeName, STR(LIB_NAME))


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to hash-table of functions with typename as the key.
//  Use when baseType requires the Targ template argument as well
#define addTemplatedToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames)                                              \
                                                                               \
    /* Add the thisType constructor function to the table */                   \
    baseType<Targ>::add##argNames##ConstructorToTable<thisType<Targ>>          \
        add##thisType##Targ##argNames##ConstructorTo##baseType##Targ##Table_;  \
                                                                               \
    bool add##thisType##Targ##argNames##To##baseType##Targ##ToC_ =             \
        addToRunTimeSelectionTableToC                                          \
        (                                                                      \
            #baseType,                                                         \
            baseType<Targ>::typeName_(),                                       \
            thisType<Targ>::typeName,                                          \
            STR(LIB_NAME)                                                      \
        )


//- Add to hash-table of functions with 'lookup' as the key.
//  Use when baseType requires the Targ template argument as well
#define addNamedTemplatedToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames,lookup)                                       \
                                                                               \
    /* Add the thisType constructor function to the table, find by lookup */   \
    baseType<Targ>::add##argNames##ConstructorToTable<thisType<Targ>>          \
        add_##lookup##_##thisType##Targ##argNames##ConstructorTo##baseType##   \
        Targ##Table_(#lookup);                                                 \
                                                                               \
    bool add_##lookup##_##thisType##Targ##argNames##To##baseType##Targ##ToC_ = \
        addToRunTimeSelectionTableToC                                          \
        (#baseType, baseType<Targ>::typeName_(), #lookup, STR(LIB_NAME))


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
