/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    RunTime selection for templated class with multiple template args.

\*---------------------------------------------------------------------------*/

#ifndef TemplatedRunTimeSelection_H
#define TemplatedRunTimeSelection_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineVariadicTemplatedRunTimeSelectionTablePtr\
(baseType,argNames,Targs...)                                                   \
    \
    /* Initialize constructorsTable pointer */                                 \
    template<>                                                                 \
    baseType<Targs >::argNames##ConstructorTable*                              \
        baseType<Targs >::argNames##ConstructorTablePtr_ = nullptr;            \
                                                                               \
    /* Specialize Construction method for the RunTime Table */                 \
    template<>                                                                 \
    void baseType<Targs >::construct##argNames##ConstructorTables()            \
    {                                                                          \
        static bool constructed = false;                                       \
        if (!constructed)                                                      \
        {                                                                      \
            constructed = true;                                                \
            baseType<Targs >::argNames##ConstructorTablePtr_                   \
                = new baseType<Targs >::argNames##ConstructorTable;            \
        }                                                                      \
    }                                                                          \
    /* Specialize Deconstruction method for the RunTime Table */               \
    template<>                                                                 \
    void baseType<Targs >::destroy##argNames##ConstructorTables()              \
    {                                                                          \
        if (baseType<Targs >::argNames##ConstructorTablePtr_)                  \
        {                                                                      \
            delete baseType<Targs >::argNames##ConstructorTablePtr_;           \
            baseType<Targs >::argNames##ConstructorTablePtr_ = nullptr;        \
        }                                                                      \
    }


#define defineVariadicTemplatedRunTimeSelectionTable(baseType,Targs...)        \
    template<>                                                                 \
    const ::Foam::word baseType<Targs >::typeName                              \
        (baseType<Targs >::typeName_());                                       \
                                                                               \
    template<>                                                                 \
    int baseType<Targs >::debug                                                \
    (                                                                          \
        ::Foam::debug::debugSwitch(baseType<Targs >::typeName_(), 0)           \
    );                                                                         \
    template<>                                                                 \
    const Foam::RegisterDebugSwitch<baseType<Targs >>                          \
    Foam::RegisterDebugSwitch<baseType<Targs >>                                \
        ::registerDebugSwitch(baseType<Targs >::typeName_());                  \
    defineVariadicTemplatedRunTimeSelectionTablePtr(baseType,dictionary,Targs);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
