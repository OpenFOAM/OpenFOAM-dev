/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Global
    Foam::globalFunctionSelectionTables

Description
    Macros to enable the easy declaration of global function selection tables.

\*---------------------------------------------------------------------------*/

#ifndef globalMemberFunctionSelectionTables_H
#define globalMemberFunctionSelectionTables_H

#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// external use:
// ~~~~~~~~~~~~~
// declare a run-time selection:
#define declareGlobalFunctionSelectionTable\
(returnType,memberFunction,argNames,argList,parList)                          \
                                                                              \
    /* Construct from argList function pointer type */                        \
    typedef returnType (*memberFunction##argNames##MemberFunctionPtr)argList; \
                                                                              \
    /* Construct from argList function table type */                          \
    typedef HashTable                                                         \
        <memberFunction##argNames##MemberFunctionPtr, word, string::hash>     \
        memberFunction##argNames##MemberFunctionTable;                        \
                                                                              \
    /* Construct from argList function pointer table pointer */               \
    extern memberFunction##argNames##MemberFunctionTable*                     \
        memberFunction##argNames##MemberFunctionTablePtr_;                    \
                                                                              \
    /* Table memberFunction called from the table add function */             \
    void construct##memberFunction##argNames##MemberFunctionTables();         \
                                                                              \
    /* Table destructor called from the table add function destructor */      \
    void destroy##memberFunction##argNames##MemberFunctionTables();           \
                                                                              \
    /* Class to add constructor from argList to table */                      \
    class add##memberFunction##argNames##GlobalMemberFunctionToTable          \
    {                                                                         \
    public:                                                                   \
                                                                              \
        add##memberFunction##argNames##GlobalMemberFunctionToTable            \
        (                                                                     \
            const word& lookup,                                               \
            memberFunction##argNames##MemberFunctionPtr function              \
        )                                                                     \
        {                                                                     \
            construct##memberFunction##argNames##MemberFunctionTables();      \
            memberFunction##argNames##MemberFunctionTablePtr_->insert         \
            (                                                                 \
                lookup,                                                       \
                function                                                      \
            );                                                                \
        }                                                                     \
                                                                              \
        ~add##memberFunction##argNames##GlobalMemberFunctionToTable()         \
        {                                                                     \
            destroy##memberFunction##argNames##MemberFunctionTables();        \
        }                                                                     \
    }


// internal use:
// constructor/destructor aid
#define defineGlobalFunctionSelectionTableConstructDestruct\
(memberFunction,argNames)                                                     \
                                                                              \
    /* Table constructor called from the table add function */                \
    void construct##memberFunction##argNames##MemberFunctionTables()          \
    {                                                                         \
        static bool constructed = false;                                      \
        if (!constructed)                                                     \
        {                                                                     \
            constructed = true;                                               \
            memberFunction##argNames##MemberFunctionTablePtr_                 \
                = new memberFunction##argNames##MemberFunctionTable;          \
        }                                                                     \
    }                                                                         \
                                                                              \
    /* Table destructor called from the table add function destructor */      \
    void destroy##memberFunction##argNames##MemberFunctionTables()            \
    {                                                                         \
        if (memberFunction##argNames##MemberFunctionTablePtr_)                \
        {                                                                     \
            delete memberFunction##argNames##MemberFunctionTablePtr_;         \
            memberFunction##argNames##MemberFunctionTablePtr_ = NULL;         \
        }                                                                     \
    }


// internal use:
// create pointer to hash-table of functions
#define defineGlobalFunctionSelectionTablePtr\
(memberFunction,argNames)                                                     \
                                                                              \
    /* Define the memberFunction table */                                     \
    memberFunction##argNames##MemberFunctionTable*                            \
        memberFunction##argNames##MemberFunctionTablePtr_ = NULL


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// external use:
// ~~~~~~~~~~~~~
// define run-time selection table
#define defineGlobalFunctionSelectionTable\
(memberFunction,argNames)                                                     \
                                                                              \
    defineGlobalFunctionSelectionTablePtr                                     \
        (memberFunction,argNames);                                            \
    defineGlobalFunctionSelectionTableConstructDestruct                       \
        (memberFunction,argNames)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
