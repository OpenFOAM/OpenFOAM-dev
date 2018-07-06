/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::RegisterSwitch

Description
    Class and registration macros for InfoSwitches and OptimisationSwitches
    to support reading from system/controlDict and dynamic update.

\*---------------------------------------------------------------------------*/

#ifndef registerSwitch_H
#define registerSwitch_H

#include "simpleRegIOobject.H"
#include "macros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RegisterSwitch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class RegisterSwitch
:
    public simpleRegIOobject
{
    Type& optSwitch_;

public:

    RegisterSwitch
    (
        void (*registryFn)(const char* name, simpleRegIOobject*),
        const char* name,
        Type& optSwitch
    )
    :
        simpleRegIOobject(registryFn, name),
        optSwitch_(optSwitch)
    {}

    virtual ~RegisterSwitch()
    {}

    virtual void readData(Foam::Istream& is)
    {
        is >> optSwitch_;
    }

    virtual void writeData(Foam::Ostream& os) const
    {
        os << optSwitch_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define registerOptSwitch(Name, Type, Switch)                                  \
    static Foam::RegisterSwitch<Type> FILE_UNIQUE(_addToOpt_)                  \
        (Foam::debug::addOptimisationObject, Name, Switch)


#define registerInfoSwitch(Name, Type, Switch)                                 \
    static Foam::RegisterSwitch<Type> FILE_UNIQUE(_addToOpt_)                  \
        (Foam::debug::addInfoObject, Name, Switch)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
