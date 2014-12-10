/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Macro definitions for debug symbols etc.

\*---------------------------------------------------------------------------*/

#ifndef debugName_H
#define debugName_H

#include "simpleRegIOobject.H"
#include "word.H"
#include "debug.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// definitions (debug information only)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//- Define the debug information, lookup as \a Name
#define registerDebugSwitchWithName(Type,Tag,Name)                            \
    class add##Tag##ToDebug                                                   \
    :                                                                         \
        public ::Foam::simpleRegIOobject                                      \
    {                                                                         \
    public:                                                                   \
        add##Tag##ToDebug(const char* name)                                   \
        :                                                                     \
            ::Foam::simpleRegIOobject(Foam::debug::addDebugObject, name)      \
        {}                                                                    \
        virtual ~add##Tag##ToDebug()                                          \
        {}                                                                    \
        virtual void readData(Foam::Istream& is)                              \
        {                                                                     \
            Type::debug = readLabel(is);                                      \
        }                                                                     \
        virtual void writeData(Foam::Ostream& os) const                       \
        {                                                                     \
            os << Type::debug;                                                \
        }                                                                     \
    };                                                                        \
    add##Tag##ToDebug add##Tag##ToDebug_(Name)


//- Register info switch (if int), lookup as \a Name
#define registerInfoSwitchWithName(Switch,Tag,Name)                           \
    class add##Tag##ToInfo                                                    \
    :                                                                         \
        public ::Foam::simpleRegIOobject                                      \
    {                                                                         \
    public:                                                                   \
        add##Tag##ToInfo(const char* name)                                    \
        :                                                                     \
            ::Foam::simpleRegIOobject(Foam::debug::addInfoObject, name)       \
        {}                                                                    \
        virtual ~add##Tag##ToInfo()                                           \
        {}                                                                    \
        virtual void readData(Foam::Istream& is)                              \
        {                                                                     \
            Switch = readLabel(is);                                           \
        }                                                                     \
        virtual void writeData(Foam::Ostream& os) const                       \
        {                                                                     \
            os << Switch;                                                     \
        }                                                                     \
    };                                                                        \
    add##Tag##ToInfo add##Tag##ToInfo_(Name)


//- Register optimisation switch (if int), lookup as \a Name
#define registerOptSwitchWithName(Switch,Tag,Name)                            \
    class add##Tag##ToOpt                                                     \
    :                                                                         \
        public ::Foam::simpleRegIOobject                                      \
    {                                                                         \
    public:                                                                   \
        add##Tag##ToOpt(const char* name)                                     \
        :                                                                     \
            ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject,name)\
        {}                                                                    \
        virtual ~add##Tag##ToOpt()                                            \
        {}                                                                    \
        virtual void readData(Foam::Istream& is)                              \
        {                                                                     \
            Switch = readLabel(is);                                           \
        }                                                                     \
        virtual void writeData(Foam::Ostream& os) const                       \
        {                                                                     \
            os << Switch;                                                     \
        }                                                                     \
    };                                                                        \
    add##Tag##ToOpt add##Tag##ToOpt_(Name)


//- Define the debug information, lookup as \a Name
#define defineDebugSwitchWithName(Type, Name, DebugSwitch)                    \
    int Type::debug(::Foam::debug::debugSwitch(Name, DebugSwitch))

//- Define the debug information
#define defineDebugSwitch(Type, DebugSwitch)                                  \
    defineDebugSwitchWithName(Type, Type::typeName_(), DebugSwitch);          \
    registerDebugSwitchWithName(Type, Type, Type::typeName_())

#ifdef __INTEL_COMPILER
//- Define the debug information for templates, lookup as \a Name
# define defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)           \
    defineDebugSwitchWithName(Type, Name, DebugSwitch)
//- Define the debug information for templates sub-classes, lookup as \a Name
# define defineTemplate2DebugSwitchWithName(Type, Name, DebugSwitch)          \
    defineDebugSwitchWithName(Type, Name, DebugSwitch)
#else
//- Define the debug information for templates, lookup as \a Name
# define defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)           \
    template<>                                                                \
    defineDebugSwitchWithName(Type, Name, DebugSwitch)
//- Define the debug information for templates sub-classes, lookup as \a Name
# define defineTemplate2DebugSwitchWithName(Type, Name, DebugSwitch)          \
    template<>                                                                \
    defineDebugSwitchWithName(Type, Name, DebugSwitch)
#endif

//- Define the debug information for templates
//  Useful with typedefs
#define defineTemplateDebugSwitch(Type, DebugSwitch)                          \
    defineTemplateDebugSwitchWithName(Type, #Type, DebugSwitch)

//- Define the debug information directly for templates
#define defineNamedTemplateDebugSwitch(Type, DebugSwitch)                     \
    defineTemplateDebugSwitchWithName(Type, Type::typeName_(), DebugSwitch)


// for templated sub-classes

//- Define the debug information for templates
//  Useful with typedefs
#define defineTemplate2DebugSwitch(Type, DebugSwitch)                         \
    defineTemplate2DebugSwitchWithName(Type, #Type, DebugSwitch)

//- Define the debug information directly for templates
#define defineNamedTemplate2DebugSwitch(Type, DebugSwitch)                    \
    defineTemplate2DebugSwitchWithName(Type, Type::typeName_(), DebugSwitch)



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
