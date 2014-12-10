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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::DataEntry<Type> > Foam::DataEntry<Type>::New
(
    const word& entryName,
    const dictionary& dict
)
{
    Istream& is(dict.lookup(entryName, false));

    token firstToken(is);

    word DataEntryType;
    if (firstToken.isWord())
    {
        // Dimensioned type default compatibility
        if (firstToken.wordToken() == entryName)
        {
            DataEntryType = "CompatibilityConstant";
        }
        else
        {
            DataEntryType = firstToken.wordToken();
        }
    }
    else
    {
//        DataEntryType = CompatibilityConstant<Type>::typeName;
        DataEntryType = "CompatibilityConstant";
    }

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(DataEntryType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("DataEntry<Type>::New(const word&, const dictionary&)")
            << "Unknown DataEntry type "
            << DataEntryType << " for DataEntry "
            << entryName << nl << nl
            << "Valid DataEntry types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<DataEntry<Type> >(cstrIter()(entryName, dict));
}


// ************************************************************************* //
