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

#include "Table.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Table<Type>::Table
(
    const word& entryName,
    const dictionary& dict
)
:
    TableBase<Type, Table<Type>>(entryName, dict)
{
    if (!dict.found(entryName))
    {
        dict.lookup("values") >> this->table_;
    }
    else
    {
        Istream& is(dict.lookup(entryName));
        word entryType(is);
        if (is.eof())
        {
            dict.lookup("values") >> this->table_;
        }
        else
        {
            is  >> this->table_;
        }
    }

    TableBase<Type, Table<Type>>::check();
}


template<class Type>
Foam::Function1s::Table<Type>::Table
(
    const word& name,
    const tableBase::boundsHandling boundsHandling,
    const word& interpolationScheme,
    const List<Tuple2<scalar, Type>>& table
)
:
    TableBase<Type, Table<Type>>
    (
        name,
        boundsHandling,
        interpolationScheme,
        table
    )
{}


template<class Type>
Foam::Function1s::Table<Type>::Table(const Table<Type>& tbl)
:
    TableBase<Type, Table<Type>>(tbl)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Table<Type>::~Table()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Table<Type>::writeEntries(Ostream& os) const
{
    TableBase<Type, Table<Type>>::writeEntries(os);

    os  << indent << "values" << this->table_
        << token::END_STATEMENT << endl;
}


// ************************************************************************* //
