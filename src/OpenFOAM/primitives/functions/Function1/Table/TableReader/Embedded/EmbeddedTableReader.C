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

#include "EmbeddedTableReader.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::TableReaders::Embedded<Type>::Embedded()
:
    TableReader<Type>()
{}


template<class Type>
Foam::TableReaders::Embedded<Type>::Embedded
(
    const word& name,
    const Function1s::unitConversions& units,
    const dictionary& dict
)
:
    TableReader<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::TableReaders::Embedded<Type>::~Embedded()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::List<Foam::Tuple2<Foam::scalar, Type>>
Foam::TableReaders::Embedded<Type>::read
(
    const Function1s::unitConversions& defaultUnits,
    const dictionary& dict,
    const word& valuesKeyword
) const
{
    Function1s::unitConversions units(defaultUnits);
    units.readIfPresent("units", dict);
    return TableReader<Type>::convertRead(units, dict.lookup(valuesKeyword));
}


template<class Type>
Foam::List<Foam::Tuple2<Foam::scalar, Type>>
Foam::TableReaders::Embedded<Type>::read
(
    const Function1s::unitConversions& units,
    Istream& is
)
{
    return TableReader<Type>::convertRead(units, is);
}


template<class Type>
void Foam::TableReaders::Embedded<Type>::write
(
    Ostream& os,
    const Function1s::unitConversions& units,
    const List<Tuple2<scalar, Type>>& table,
    const word& valuesKeyword
) const
{
    writeEntry
    (
        os,
        valuesKeyword,
        TableReader<Type>::convertWrite(units, table)
    );
}


// ************************************************************************* //
