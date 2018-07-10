/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Class
    Foam::memberFunctionSelectionTables

Description
    Macros to enable the easy declaration of member function selection tables.

\*---------------------------------------------------------------------------*/

#include "token.H"

#ifndef memberFunctionSelectionTables_H
#define memberFunctionSelectionTables_H

#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Declare a run-time selection:
#define declareMemberFunctionSelectionTable(                                   \
    returnType,baseType,memberFunction,argNames,argList,parList)               \
                                                                               \
    /* Construct from argList function pointer type */                         \
    typedef returnType (*memberFunction##argNames##MemberFunctionPtr)argList;  \
                                                                               \
    /* Construct from argList function table type */                           \
    typedef HashTable                                                          \
        <memberFunction##argNames##MemberFunctionPtr, word, string::hash>      \
        memberFunction##argNames##MemberFunctionTable;                         \
                                                                               \
    /* Construct from argList function pointer table pointer */                \
    static memberFunction##argNames##MemberFunctionTable*                      \
        memberFunction##argNames##MemberFunctionTablePtr_;                     \
                                                                               \
    /* Class to add constructor from argList to table */                       \
    template<class baseType##Type>                                             \
    class add##memberFunction##argNames##MemberFunctionToTable                 \
    {                                                                          \
    public:                                                                    \
                                                                               \
        add##memberFunction##argNames##MemberFunctionToTable                   \
        (                                                                      \
            const word& lookup = baseType##Type::typeName                      \
        )                                                                      \
        {                                                                      \
            construct##memberFunction##argNames##MemberFunctionTables();       \
            memberFunction##argNames##MemberFunctionTablePtr_->insert          \
            (                                                                  \
                lookup,                                                        \
                baseType##Type::memberFunction                                 \
            );                                                                 \
        }                                                                      \
                                                                               \
        ~add##memberFunction##argNames##MemberFunctionToTable()                \
        {                                                                      \
            destroy##memberFunction##argNames##MemberFunctionTables();         \
        }                                                                      \
    };                                                                         \
                                                                               \
    /* Table memberFunction called from the table add function */              \
    static void construct##memberFunction##argNames##MemberFunctionTables();   \
                                                                               \
    /* Table destructor called from the table add function destructor */       \
    static void destroy##memberFunction##argNames##MemberFunctionTables()


// Constructor aid
#define defineMemberFunctionSelectionTableMemberFunction(                      \
    baseType,memberFunction,argNames)                                          \
                                                                               \
    /* Table memberFunction called from the table add function */              \
    void baseType::construct##memberFunction##argNames##MemberFunctionTables() \
    {                                                                          \
        static bool constructed = false;                                       \
        if (!constructed)                                                      \
        {                                                                      \
            constructed = true;                                                \
            baseType::memberFunction##argNames##MemberFunctionTablePtr_        \
                = new baseType::memberFunction##argNames##MemberFunctionTable; \
        }                                                                      \
    }


// Destructor aid
#define defineMemberFunctionSelectionTableDestructor(                          \
    baseType,memberFunction,argNames)                                          \
                                                                               \
    /* Table destructor called from the table add function destructor */       \
    void baseType::destroy##memberFunction##argNames##MemberFunctionTables()   \
    {                                                                          \
        if (baseType::memberFunction##argNames##MemberFunctionTablePtr_)       \
        {                                                                      \
            delete baseType::memberFunction##argNames##MemberFunctionTablePtr_;\
            baseType::memberFunction##argNames##MemberFunctionTablePtr_ =      \
                nullptr;                                                       \
        }                                                                      \
    }


// Create pointer to hash-table of functions
#define defineMemberFunctionSelectionTablePtr(baseType,memberFunction,argNames)\
                                                                               \
    /* Define the memberFunction table */                                      \
    baseType::memberFunction##argNames##MemberFunctionTable*                   \
        baseType::memberFunction##argNames##MemberFunctionTablePtr_ = nullptr



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define run-time selection table
#define defineMemberFunctionSelectionTable(baseType,memberFunction,argNames)   \
                                                                               \
    defineMemberFunctionSelectionTablePtr                                      \
        (baseType,memberFunction,argNames);                                    \
    defineMemberFunctionSelectionTableMemberFunction                           \
        (baseType,memberFunction,argNames)                                     \
    defineMemberFunctionSelectionTableDestructor                               \
        (baseType,memberFunction,argNames)


//- Define run-time selection table for template classes
//  use when baseType doesn't need a template argument (eg, is a typedef)
#define defineTemplateMemberFunctionSelectionTable(                            \
    baseType,memberFunction,argNames)                                          \
                                                                               \
    template<>                                                                 \
    defineMemberFunctionSelectionTablePtr                                      \
        (baseType,memberFunction,argNames);                                    \
    template<>                                                                 \
    defineMemberFunctionSelectionTableMemberFunction                           \
        (baseType,memberFunction,argNames)                                     \
    template<>                                                                 \
    defineMemberFunctionSelectionTableDestructor                               \
        (baseType,memberFunction,argNames)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Constructor aid: use when baseType requires the Targ template argument
#define defineTemplatedMemberFunctionSelectionTableMemberFunction(             \
    baseType,memberFunction,argNames,Targ)                                     \
                                                                               \
    /* Table memberFunction called from the table add function */              \
    void baseType<Targ>::construct##memberFunction##argNames##                 \
    MemberFunctionTables()                                                     \
    {                                                                          \
        static bool constructed = false;                                       \
        if (!constructed)                                                      \
        {                                                                      \
            constructed = true;                                                \
            baseType<Targ>::memberFunction##argNames##MemberFunctionTablePtr_  \
                = new baseType<Targ>::memberFunction##argNames##               \
                  MemberFunctionTable;                                         \
        }                                                                      \
    }


// Destructor aid
// use when baseType requires the Targ template argument
#define defineTemplatedMemberFunctionSelectionTableDestructor(                 \
    baseType,memberFunction,argNames,Targ)                                     \
                                                                               \
    /* Table destructor called from the table add function destructor */       \
    void baseType<Targ>::destroy##memberFunction##argNames##                   \
    MemberFunctionTables()                                                     \
    {                                                                          \
        if                                                                     \
        (                                                                      \
            baseType<Targ>::memberFunction##argNames##MemberFunctionTablePtr_  \
        )                                                                      \
        {                                                                      \
            delete baseType<Targ>::memberFunction##argNames##                  \
                MemberFunctionTablePtr_;                                       \
            baseType<Targ>::memberFunction##argNames##                         \
                MemberFunctionTablePtr_ = nullptr;                             \
        }                                                                      \
    }


// Create pointer to hash-table of functions
// use when baseType requires the Targ template argument
#define defineTemplatedMemberFunctionSelectionTablePtr(                        \
    baseType,memberFunction,argNames,Targ)                                     \
                                                                               \
    /* Define the memberFunction table */                                      \
    baseType<Targ>::memberFunction##argNames##MemberFunctionTable*             \
        baseType<Targ>::memberFunction##argNames##MemberFunctionTablePtr_ =    \
        nullptr


//- Define run-time selection table for template classes
//  use when baseType requires the Targ template argument
#define defineTemplatedMemberFunctionSelectionTable(                           \
    baseType,memberFunction,argNames,Targ)                                     \
                                                                               \
    template<>                                                                 \
    defineTemplatedMemberFunctionSelectionTablePtr                             \
        (baseType,memberFunction,argNames,Targ);                               \
    template<>                                                                 \
    defineTemplatedMemberFunctionSelectionTableMemberFunction                  \
        (baseType,memberFunction,argNames,Targ)                                \
    template<>                                                                 \
    defineTemplatedMemberFunctionSelectionTableDestructor                      \
        (baseType,memberFunction,argNames,Targ)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
