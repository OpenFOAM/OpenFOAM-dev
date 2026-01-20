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

#include "TableReader.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Coordinate, class Value>
void Foam::TableReader<Coordinate, Value>::convertRead
(
    const Function1s::unitConversions& units,
    List<Tuple2<Coordinate, Value>>& table
) const
{
    forAll(table, i)
    {
        table[i].first() = units.x.toStandard(table[i].first());
        table[i].second() = units.value.toStandard(table[i].second());
    }
}


template<class Coordinate, class Value>
Foam::List<Foam::Tuple2<Coordinate, Value>>
Foam::TableReader<Coordinate, Value>::convertRead
(
    const Function1s::unitConversions& units,
    const List<Tuple2<Coordinate, Value>>& table
) const
{
    List<Tuple2<Coordinate, Value>> tableCopy(table);
    convertRead(units, tableCopy);
    return tableCopy;
}


template<class Coordinate, class Value>
Foam::List<Foam::Tuple2<Coordinate, Value>>
Foam::TableReader<Coordinate, Value>::convertWrite
(
    const Function1s::unitConversions& units,
    const List<Tuple2<Coordinate, Value>>& table
) const
{
    List<Tuple2<Coordinate, Value>> tableCopy(table);

    forAll(tableCopy, i)
    {
        tableCopy[i].first() = units.x.toUser(tableCopy[i].first());
        tableCopy[i].second() = units.value.toUser(tableCopy[i].second());
    }

    return tableCopy;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Coordinate, class Value>
Foam::TableReader<Coordinate, Value>::TableReader()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Coordinate, class Value>
Foam::TableReader<Coordinate, Value>::~TableReader()
{}


// ************************************************************************* //
