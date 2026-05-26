/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "Constant.H"
#include "unitSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Constant<Type>::Constant
(
    const word& name,
    const Type& val
)
:
    FieldFunction1<Type, Constant<Type>>(name),
    value_(val)
{}


template<class Type>
Foam::Function1s::Constant<Type>::Constant
(
    const word& name,
    const unitSets& units,
    const dictionary& dict
)
:
    FieldFunction1<Type, Constant<Type>>(name),
    value_(dict.lookup<Type>("value", typeUnits<Type>(units.value)))
{}


template<class Type>
Foam::Function1s::Constant<Type>::Constant
(
    const word& name,
    const unitSets& units,
    Istream& is
)
:
    FieldFunction1<Type, Constant<Type>>(name),
    value_(readAndConvert<Type>(is, typeUnits<Type>(units.value)))
{}


template<class Type>
Foam::Function1s::Constant<Type>::Constant(const Constant<Type>& cnst)
:
    FieldFunction1<Type, Constant<Type>>(cnst),
    value_(cnst.value_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Constant<Type>::~Constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Constant<Type>::write
(
    Ostream& os,
    const unitSets& units
) const
{
    writeEntry(os, "value", typeUnits<Type>(units.value), value_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Constant<Type>::operator=(const Constant<Type>& cnst)
{
    FieldFunction1<Type, Constant<Type>>::operator=(cnst);

    value_ = cnst.value_;
}


// ************************************************************************* i/
