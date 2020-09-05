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
    Macros to add child classes to the parent's RunTime selection table.

\*---------------------------------------------------------------------------*/

#ifndef addToTemplatedRunTimeSelection_H
#define addToTemplatedRunTimeSelection_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeTemplatedModel(baseType,childType,TargName,ns,Targs...)            \
                                                                               \
    /* Define child type name */                                               \
    template<>                                                                 \
    const ::Foam::word ns::childType<Targs >::                                 \
        typeName(ns::childType<Targs >::typeName_());                          \
                                                                               \
    /* Define child debug switch */                                            \
    template<>                                                                 \
    int ns::childType<Targs >::debug                                           \
    (                                                                          \
        ::Foam::debug::debugSwitch(ns::childType<Targs >::typeName_(), 0)      \
    );                                                                         \
    /* Register the debug switch */                                            \
    template<>                                                                 \
    const Foam::RegisterDebugSwitch<ns::childType<Targs >>                     \
    Foam::RegisterDebugSwitch<ns::childType<Targs >>                           \
        ::registerDebugSwitch(ns::childType<Targs >::typeName_());             \
                                                                               \
    /* Add child to the base's run-time selection tabe */                      \
    namespace ns {                                                             \
    baseType<Targs >::adddictionaryConstructorToTable<childType<Targs >>       \
    add##phaseType##TargName##argNames##ConstructorTo##baseType##Targ##Table_; \
    }
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
