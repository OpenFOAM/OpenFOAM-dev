/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "Function12.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Function12<Type>::Function12
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction2<Type, Function12<Type>>(name),
    index_(),
    f_()
{
    static const word name1 = "value1", name2 = "value2";
    const bool found1 = dict.found(name1), found2 = dict.found(name2);

    if (found1 && found2)
    {
        FatalIOErrorInFunction(dict)
            << "Both keywords " << name1 << " and " << name2
            << " are defined in dictionary " << dict.name()
            << exit(FatalError);
    }

    if (!found1 && !found2)
    {
        FatalIOErrorInFunction(dict)
            << "Neither keywords " << name1 << " or " << name2
            << " are defined in dictionary " << dict.name()
            << exit(FatalError);
    }

    index_ = found2;

    f_ =
        Function1<Type>::New
        (
            found2 ? name2 : name1,
            {found2 ? units.y : units.x, units.value},
            dict
        );
}


template<class Type>
Foam::Function2s::Function12<Type>::Function12(const Function12<Type>& se)
:
    FieldFunction2<Type, Function12<Type>>(se),
    index_(se.index_),
    f_(se.f_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Function12<Type>::~Function12()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function2s::Function12<Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, {index_ ? units.y : units.x, units.value}, f_());
}


// ************************************************************************* //
