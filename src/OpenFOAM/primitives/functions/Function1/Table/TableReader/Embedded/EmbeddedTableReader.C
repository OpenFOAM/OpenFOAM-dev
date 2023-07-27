/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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
Foam::TableReaders::Embedded<Type>::read(const dictionary& dict) const
{
    return dict.lookup<List<Tuple2<scalar, Type>>>("values");
}


template<class Type>
void Foam::TableReaders::Embedded<Type>::write
(
    Ostream& os,
    const List<Tuple2<scalar, Type>>& table
) const
{
    writeEntry(os, "values", table);
}


// ************************************************************************* //
