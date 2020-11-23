/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
void Foam::TableFileReader<Type>::read
(
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
            << "Cannot open file" << fName_ << nl
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::TableFileReader<Type>::TableFileReader
(
    const dictionary& dict
)
:
    TableReader<Type>(dict),
    fName_(dict.lookup("file"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::TableFileReader<Type>::~TableFileReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::TableFileReader<Type>::write
(
    Ostream& os,
    const List<Tuple2<scalar, Type>>& table
) const
{
    writeEntry(os, "format", this->type());
    writeEntry(os, "file", fName_);
}


// ************************************************************************* //
