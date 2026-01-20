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

#include "FoamTableReader.H"
#include "tokenTupleN.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Coordinate, class Value>
void Foam::TableReaders::Foam<Coordinate, Value>::read
(
    ISstream& is,
    List<Tuple2<Coordinate, Value>>& data
) const
{
    List<tokenTupleN> dataStr(is);

    data.resize(dataStr.size());

    for (label i = 0; i < dataStr.size(); ++ i)
    {
        data[i].first() = dataStr[i].get<Coordinate>(is, columns_.first());
        data[i].second() = dataStr[i].get<Value>(is, columns_.second());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Coordinate, class Value>
Foam::TableReaders::Foam<Coordinate, Value>::Foam
(
    const word& name,
    const Function1s::unitConversions& units,
    const dictionary& dict
)
:
    TableFileReader<Coordinate, Value>(units, dict),
    columns_(dict.lookupOrDefault<labelPair>("columns", labelPair(0, 1)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Coordinate, class Value>
Foam::TableReaders::Foam<Coordinate, Value>::~Foam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Coordinate, class Value>
void Foam::TableReaders::Foam<Coordinate, Value>::write
(
    Ostream& os,
    const Function1s::unitConversions& units,
    const List<Tuple2<Coordinate, Value>>& table,
    const word&
) const
{
    TableFileReader<Coordinate, Value>::write(os, units, table);

    writeEntry(os, "columns", columns_);
}


// ************************************************************************* //
