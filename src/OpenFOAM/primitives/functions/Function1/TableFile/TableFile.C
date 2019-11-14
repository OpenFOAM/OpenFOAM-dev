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

#include "TableFile.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::TableFile<Type>::TableFile
(
    const word& entryName,
    const dictionary& dict
)
:
    TableBase<Type, TableFile<Type>>(entryName, dict),
    fName_(dict.lookup("file")),
    reader_
    (
        TableReader<Type>::New
        (
            dict.lookupOrDefault<word>
            (
                "format",
                TableReaders::Foam<Type>::typeName
            ),
            dict
        )
    )
{
    reader_()(fName_, this->table_);

    TableBase<Type, TableFile<Type>>::check();
}


template<class Type>
Foam::Function1s::TableFile<Type>::TableFile(const TableFile<Type>& tbl)
:
    TableBase<Type, TableFile<Type>>(tbl),
    fName_(tbl.fName_),
    reader_(tbl.reader_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::TableFile<Type>::~TableFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::TableFile<Type>::writeEntries(Ostream& os) const
{
    TableBase<Type, TableFile<Type>>::writeEntries(os);

    writeEntry(os, "file", fName_);
    writeEntry(os, "format", reader_->type());
    reader_->write(os);
}


// ************************************************************************* //
