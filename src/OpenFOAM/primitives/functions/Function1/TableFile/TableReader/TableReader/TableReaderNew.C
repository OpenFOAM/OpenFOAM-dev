/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
#include "FoamTableReader.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::TableReader<Type>> Foam::TableReader<Type>::New
(
    const word& name,
    const dictionary& dict,
    List<Tuple2<scalar, Type>>& table
)
{
    if (dict.found("format"))
    {
        const word readerType(dict.lookup("format"));

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

        return autoPtr<TableReader<Type>>(cstrIter()(name, dict, table));
    }
    else
    {
        if (dict.found("file"))
        {
            return autoPtr<TableReader<Type>>
            (
                new TableReaders::Foam<Type>(name, dict, table)
            );
        }
        else
        {
            return autoPtr<TableReader<Type>>
            (
                new TableReaders::Embedded<Type>(name, dict, table)
            );
        }
    }
}


// ************************************************************************* //
