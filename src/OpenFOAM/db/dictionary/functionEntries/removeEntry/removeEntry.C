/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "removeEntry.H"
#include "dictionary.H"
#include "stringListOps.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionEntries::removeEntry::typeName
(
    Foam::functionEntries::removeEntry::typeName_()
);

// Don't lookup the debug switch here as the debug switch dictionary
// might include removeEntry
int Foam::functionEntries::removeEntry::debug(0);

namespace Foam
{
namespace functionEntries
{
    addToMemberFunctionSelectionTable
    (
        functionEntry,
        removeEntry,
        execute,
        dictionaryIstream
    );
}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::removeEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    wordList   dictKeys = parentDict.toc();
    wordReList patterns = readList<wordRe>(is);

    labelList indices = findStrings(patterns, dictKeys);

    forAll(indices, indexI)
    {
        parentDict.remove(dictKeys[indices[indexI]]);
    }

    return true;
}

// ************************************************************************* //
