/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
#include "fileOperation.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::TableReader<Type>> Foam::TableReader<Type>::New
(
    const word& readerType,
    const dictionary& dict
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(readerType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown reader type " << readerType
            << nl << nl
            << "Valid reader types : " << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<TableReader<Type>>(cstrIter()(dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::TableReader<Type>::TableReader(const dictionary&)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::TableReader<Type>::~TableReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::TableReader<Type>::write(Ostream& os) const
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
template<class TableType>
void Foam::TableReader<Type>::operator()
(
    const fileName& fName,
    TableType& data
) const
{
    // Expand the file
    fileName fNameExpanded(fName);
    fNameExpanded.expand();

    // Open a stream and check it
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(fNameExpanded));
    ISstream& is = isPtr();
    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot open file" << fName << nl
            << exit(FatalIOError);
    }

    // Read data from the stream
    read(is, data);

    // Check something was read
    if (data.empty())
    {
        FatalIOErrorInFunction(is)
            << "Table read from " << fName << " is empty" << nl
            << exit(FatalIOError);
    }
}


// ************************************************************************* //
