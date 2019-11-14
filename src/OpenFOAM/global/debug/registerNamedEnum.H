/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::RegisterNamedEnum

Description
    Class and registration macros for info and optimisation named enumerations
    to support reading from system/controlDict and dynamic update.

\*---------------------------------------------------------------------------*/

#ifndef registerNamedEnum_H
#define registerNamedEnum_H

#include "simpleRegIOobject.H"
#include "macros.H"
#include "NamedEnum.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class RegisterNamedEnum Declaration
\*---------------------------------------------------------------------------*/

template<class NamedEnumType>
class RegisterNamedEnum;

template<class Enum, unsigned int nEnum>
class RegisterNamedEnum<NamedEnum<Enum, nEnum>>
:
    public simpleRegIOobject
{
    const NamedEnum<Enum, nEnum>& enumNames_;

    Enum& enumRef_;

public:

    RegisterNamedEnum
    (
        void (*registryFn)(const char* name, simpleRegIOobject*),
        const word& name,
        const NamedEnum<Enum, nEnum>& enumNames,
        Enum& enumRef
    )
    :
        ::Foam::simpleRegIOobject(registryFn, name.c_str()),
        enumNames_(enumNames),
        enumRef_(enumRef)
    {}

    virtual ~RegisterNamedEnum()
    {}

    virtual void readData(Istream& is)
    {
        enumRef_ = enumNames_.read(is);
    }

    virtual void writeData(Ostream& os) const
    {
        os  << enumNames_[enumRef_];
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define registerOptNamedEnum(Name, EnumNames, EnumRef)                         \
    static Foam::RegisterNamedEnum<std::remove_cv<decltype(EnumNames)>::type>  \
        FILE_UNIQUE(_addToOpt_)                                                \
        (Foam::debug::addOptimisationObject, Name, EnumNames, EnumRef)

#define registerInfoNamedEnum(Name, EnumNames, EnumRef)                        \
    static Foam::RegisterNamedEnum<std::remove_cv<decltype(EnumNames)>::type>  \
        FILE_UNIQUE(_addToOpt_)                                                \
        (Foam::debug::addInfoObject, Name, EnumNames, EnumRef)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
