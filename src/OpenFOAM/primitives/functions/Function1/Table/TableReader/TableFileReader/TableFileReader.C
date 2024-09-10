/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "TableFileReader.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::Function1s::unitConversions>
Foam::TableFileReader<Type>::readUnits
(
    const Function1s::unitConversions& defaultUnits,
    const dictionary& dict
) const
{
    if (dict.found("units"))
    {
        autoPtr<Function1s::unitConversions> unitsPtr
        (
            new Function1s::unitConversions(defaultUnits)
        );
        unitsPtr->readIfPresent("units", dict);
        return unitsPtr;
    }
    else
    {
        return autoPtr<Function1s::unitConversions>(nullptr);
    }
}


template<class Type>
void Foam::TableFileReader<Type>::read
(
    const Function1s::unitConversions& defaultUnits,
    const dictionary& dict,
    List<Tuple2<scalar, Type>>& table
) const
{
    // Expand the file
    fileName fNameExpanded(fName_);
    fNameExpanded.expand();

    // Open a stream and check it
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(fNameExpanded));
    ISstream& is = isPtr();
    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot open file " << fName_ << nl
            << exit(FatalIOError);
    }

    // Read data from the stream
    read(is, table);

    // Check something was read
    if (table.empty())
    {
        FatalIOErrorInFunction(is)
            << "Table read from " << fName_ << " is empty" << nl
            << exit(FatalIOError);
    }

    // Convert units
    TableReader<Type>::convertRead
    (
        unitsPtr_.valid() ? unitsPtr_() : defaultUnits,
        table
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::TableFileReader<Type>::TableFileReader
(
    const Function1s::unitConversions& defaultUnits,
    const dictionary& dict
)
:
    TableReader<Type>(),
    fName_(dict.lookup("file")),
    unitsPtr_(readUnits(defaultUnits, dict))
{}


template<class Type>
Foam::TableFileReader<Type>::TableFileReader
(
    const TableFileReader<Type>& tfr
)
:
    TableReader<Type>(tfr),
    fName_(tfr.fName_),
    unitsPtr_(tfr.unitsPtr_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::TableFileReader<Type>::~TableFileReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::List<Foam::Tuple2<Foam::scalar, Type>>
Foam::TableFileReader<Type>::read
(
    const Function1s::unitConversions& units,
    const dictionary& dict,
    const word&
) const
{
    List<Tuple2<scalar, Type>> data;
    read(units, dict, data);
    return data;
}


template<class Type>
void Foam::TableFileReader<Type>::write
(
    Ostream& os,
    const Function1s::unitConversions& units,
    const List<Tuple2<scalar, Type>>& table,
    const word&
) const
{
    writeEntry(os, "format", this->type());
    writeEntry(os, "file", fName_);

    if (unitsPtr_.valid())
    {
        writeEntry(os, "units", unitsPtr_());
    }
}


// ************************************************************************* //
