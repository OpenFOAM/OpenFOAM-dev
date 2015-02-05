/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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
    Macro definitions for info switches.

\*---------------------------------------------------------------------------*/

#ifndef registerInfoSwitch_H
#define registerInfoSwitch_H

#include "simpleRegIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
