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

#include "DataEntry.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const TableFile<Type>& tbl
)
{
    os  << static_cast<const DataEntry<Type>&>(tbl)
        << static_cast<const TableBase<Type>&>(tbl);

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const TableFile<Type>&)"
    );

    return os;
}


template<class Type>
void Foam::TableFile<Type>::writeData(Ostream& os) const
{
    DataEntry<Type>::writeData(os);

    os  << token::END_STATEMENT << nl
        << indent << word(type() + "Coeffs") << nl
        << indent << token::BEGIN_BLOCK << nl << incrIndent;

    // Note: for TableBase write the dictionary entries it needs but not
    // the values themselves
    TableBase<Type>::writeEntries(os);

    os.writeKeyword("fileName")<< fName_ << token::END_STATEMENT << nl;
    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// ************************************************************************* //
