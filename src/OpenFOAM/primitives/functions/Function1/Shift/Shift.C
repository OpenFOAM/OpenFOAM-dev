/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "Shift.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Shift<Type>::Shift
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction1<Type, Shift<Type>>(name),
    start_(dict.lookup<scalar>("start", units.x)),
    value_(Function1<Type>::New("value", units, dict))
{}


template<class Type>
Foam::Function1s::Shift<Type>::Shift(const Shift<Type>& se)
:
    FieldFunction1<Type, Shift<Type>>(se),
    start_(se.start_),
    value_(se.value_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Shift<Type>::~Shift()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Shift<Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, "start", units.x, start_);
    writeEntry(os, units, value_());
}


// ************************************************************************* //
