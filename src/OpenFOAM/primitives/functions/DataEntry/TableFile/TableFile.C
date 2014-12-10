/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
Foam::TableFile<Type>::TableFile(const word& entryName, const dictionary& dict)
:
    DataEntry<Type>(entryName),
    TableBase<Type>(entryName, dict.subDict(type() + "Coeffs")),
    fName_("none")
{
    const dictionary coeffs(dict.subDict(type() + "Coeffs"));
    coeffs.lookup("fileName") >> fName_;

    if (coeffs.found("dimensions"))
    {
        coeffs.lookup("dimensions") >> this->dimensions_;
    }

    fileName expandedFile(fName_);
    IFstream is(expandedFile.expand());

    if (!is.good())
    {
        FatalIOErrorIn
        (
            "TableFile<Type>::TableFile(const word&, const dictionary&)",
            is
        )   << "Cannot open file." << exit(FatalIOError);
    }

    is  >> this->table_;

    TableBase<Type>::check();
}


template<class Type>
Foam::TableFile<Type>::TableFile(const TableFile<Type>& tbl)
:
    DataEntry<Type>(tbl),
    TableBase<Type>(tbl),
    fName_(tbl.fName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::TableFile<Type>::~TableFile()
{}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "TableFileIO.C"


// ************************************************************************* //
