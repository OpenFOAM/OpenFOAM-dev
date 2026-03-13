/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "DimensionedFieldFunction.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class DimensionedFieldType>
Foam::DimensionedFieldFunction<DimensionedFieldType>::DimensionedFieldFunction
(
    const dictionary& dict,
    DimensionedFieldType& field
)
:
    field_(field)
{}


template<class DimensionedFieldType>
Foam::DimensionedFieldFunction<DimensionedFieldType>::DimensionedFieldFunction
(
    const DimensionedFieldFunction& dff,
    DimensionedFieldType& field
)
:
    field_(field)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class DimensionedFieldType>
Foam::autoPtr<Foam::DimensionedFieldFunction<DimensionedFieldType>>
Foam::DimensionedFieldFunction<DimensionedFieldType>::New
(
    const dictionary& dict,
    DimensionedFieldType& field
)
{
    const word type(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown DimensionedFieldFunction type " << type
            << " for field " << field.name() << nl << nl
            << "Valid DimensionedFieldFunction types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return cstrIter()(dict, field);
}


// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

template<class DimensionedFieldType>
void Foam::writeEntry
(
    Ostream& os,
    const word& key,
    const DimensionedFieldFunction<DimensionedFieldType>& f
)
{
    writeKeyword(os, key)
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

    writeEntry(os, "type", f.type());
    f.write(os);

    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// ************************************************************************* //
