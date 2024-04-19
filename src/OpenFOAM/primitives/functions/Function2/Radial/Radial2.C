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

#include "Radial2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Radial<Type>::Radial
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction2<Type, Radial<Type>>(name),
    value_
    (
        Function1<Type>::New("value", {units.x + units.y, units.value}, dict)
    )
{}


template<class Type>
Foam::Function2s::Radial<Type>::Radial(const Radial<Type>& se)
:
    FieldFunction2<Type, Radial<Type>>(se),
    value_(se.value_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Radial<Type>::~Radial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function2s::Radial<Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, {units.x + units.y, units.value}, value_());
}


// ************************************************************************* //
