/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

template<class Type>
Foam::demandDrivenEntry<Type>::demandDrivenEntry
(
    const dictionary& dict,
    const Type& value
)
:
    dict_(dict),
    keyword_("unknown-keyword"),
    value_(value),
    stored_(true)
{}


template<class Type>
Foam::demandDrivenEntry<Type>::demandDrivenEntry
(
    const dictionary& dict,
    const word& keyword
)
:
    dict_(dict),
    keyword_(keyword),
    value_(Zero),
    stored_(false)
{}


template<class Type>
Foam::demandDrivenEntry<Type>::demandDrivenEntry
(
    const dictionary& dict,
    const word& keyword,
    const Type& defaultValue,
    const bool readIfPresent
)
:
    dict_(dict),
    keyword_(keyword),
    value_(defaultValue),
    stored_(true)
{
    if (readIfPresent)
    {
        dict_.readIfPresent<Type>(keyword, value_);
    }
}


template<class Type>
Foam::demandDrivenEntry<Type>::demandDrivenEntry(const demandDrivenEntry& dde)
:
    dict_(dde.dict_),
    keyword_(dde.keyword_),
    value_(dde.value_),
    stored_(dde.stored_)
{}


// ************************************************************************* //
