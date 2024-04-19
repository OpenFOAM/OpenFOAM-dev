/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2024 OpenFOAM Foundation
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

#include "Constant2.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Type Foam::Function2s::Constant<Type>::readValue
(
    const unitConversions& defaultUnits,
    Istream& is
)
{
    // Read the units if they are before the value
    unitConversion units(defaultUnits.value);
    const bool haveUnits = units.readIfPresent(is);

    // Read the value
    const Type value = pTraits<Type>(is);

    // Read the units if they are after the value
    if (!haveUnits && !is.eof())
    {
        units.readIfPresent(is);
    }

    // Modify the value by the unit conversion and return
    return units.toStandard(value);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Constant<Type>::Constant
(
    const word& name,
    const Type& val
)
:
    FieldFunction2<Type, Constant<Type>>(name),
    value_(val)
{}


template<class Type>
Foam::Function2s::Constant<Type>::Constant
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction2<Type, Constant<Type>>(name),
    value_(dict.lookup<Type>("value", units.value))
{}


template<class Type>
Foam::Function2s::Constant<Type>::Constant
(
    const word& name,
    const unitConversions& units,
    Istream& is
)
:
    FieldFunction2<Type, Constant<Type>>(name),
    value_(readValue(units, is))
{}


template<class Type>
Foam::Function2s::Constant<Type>::Constant(const Constant<Type>& cnst)
:
    FieldFunction2<Type, Constant<Type>>(cnst),
    value_(cnst.value_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Constant<Type>::~Constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function2s::Constant<Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, "value", units.value, value_);
}


// ************************************************************************* i/
