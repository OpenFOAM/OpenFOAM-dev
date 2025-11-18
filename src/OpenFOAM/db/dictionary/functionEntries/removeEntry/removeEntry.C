/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
#include "stringListOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(removeEntry, 0);
    addToRunTimeSelectionTable(functionEntry, removeEntry, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tokenList Foam::functionEntries::removeEntry::readArgOrList
(
    const functionName& functionType,
    Istream& is
) const
{
    tokenList argList;

    // Read the next token to check for '('
    // in case the optional arguments start on the next line
    token currToken(is);

    if
    (
        currToken.isPunctuation()
     && currToken.pToken() == token::BEGIN_LIST
    )
    {
        do
        {
            argList.append(currToken);
        }
        while
        (
            !(currToken == token::END_LIST)
         && !is.read(currToken).bad()
         && currToken.good()
        );

        if
        (
            !currToken.isPunctuation()
         || currToken.pToken() != token::END_LIST
        )
        {
            FatalIOErrorInFunction(is)
                << "Unclosed argument list " << currToken
                << " in functionEntry " << functionType
                << exit(FatalIOError);
        }
    }
    else
    {
        argList.append(currToken);
    }

    return argList;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::removeEntry::removeEntry
(
    const dictionary& parentDict,
    Istream& is
)
:
    functionEntry(typeName, parentDict, is, readArgOrList(typeName, is)),
    patterns_(readList<wordRe>(stream()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::removeEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    const wordList dictKeys = contextDict.toc();
    const labelList indices = findStrings(patterns_, dictKeys);

    forAll(indices, i)
    {
        contextDict.remove(dictKeys[indices[i]]);
    }

    return true;
}


// ************************************************************************* //
