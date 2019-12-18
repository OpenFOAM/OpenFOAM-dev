/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "negEntry.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(negEntry, 0);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        negEntry,
        execute,
        dictionaryIstream
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        negEntry,
        execute,
        primitiveEntryIstream
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::string Foam::functionEntries::negEntry::negateVariable
(
    const dictionary& parentDict,
    Istream& is
)
{
    // Read variable name as a word including the '$'
    const variable var(is);

    // Strip the leading '$' from the variable name
    const string varName = var(1, var.size() - 1);

    // Lookup the variable name in the parent dictionary....
    const entry* ePtr = parentDict.lookupScopedEntryPtr(varName, true, false);

    if (ePtr && ePtr->isStream())
    {
        const token variable(ePtr->stream());

        // Convert to a string
        OStringStream os(is.format());
        os << variable;
        const string str(os.str());

        // Negate
        if (str[0] == '-')
        {
            return str(1, str.size() - 1);
        }
        else
        {
            return '-' + str;
        }
    }
    else
    {
        FatalIOErrorInFunction
        (
            parentDict
        )   << "Illegal dictionary variable name " << varName << endl
            << "Valid dictionary entries are " << parentDict.toc()
            << exit(FatalIOError);

        return string::null;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::negEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    return insert(parentDict, negateVariable(parentDict, is));
}


bool Foam::functionEntries::negEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& thisEntry,
    Istream& is
)
{
    return insert(parentDict, thisEntry, negateVariable(parentDict, is));
}


// ************************************************************************* //
