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

#include "inputModeEntry.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::functionEntries::inputModeEntry::inputMode
    Foam::functionEntries::inputModeEntry::mode_(MERGE);

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(inputModeEntry, 0);
    addToRunTimeSelectionTable(functionEntry, inputModeEntry, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::inputModeEntry::inputModeEntry
(
    const label lineNumber,
    const dictionary& parentDict,
    Istream& is
)
:
    functionEntry(typeName, lineNumber, parentDict, is, token(is))
{
    if (!operator[](0).isWord())
    {
        FatalIOErrorInFunction(is)
            << "Expected a word, found " << operator[](0)
            << " while reading function " << typeName
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionEntries::inputModeEntry::clear()
{
    mode_ = MERGE;
}


bool Foam::functionEntries::inputModeEntry::merge()
{
    return mode_ == MERGE;
}


bool Foam::functionEntries::inputModeEntry::overwrite()
{
    return mode_ == OVERWRITE;
}


bool Foam::functionEntries::inputModeEntry::protect()
{
    return mode_ == PROTECT;
}

bool Foam::functionEntries::inputModeEntry::error()
{
    return mode_ == ERROR;
}


bool Foam::functionEntries::inputModeEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    clear();

    const word& mode = operator[](0).wordToken();

    if (mode == "merge" || mode == "default")
    {
        mode_ = MERGE;
    }
    else if (mode == "overwrite")
    {
        mode_ = OVERWRITE;
    }
    else if (mode == "protect")
    {
        mode_ = PROTECT;
    }
    else if (mode == "warn")
    {
        mode_ = WARN;
    }
    else if (mode == "error")
    {
        mode_ = ERROR;
    }
    else
    {
        WarningInFunction
            << "unsupported input mode '" << mode
            << "' ... defaulting to 'merge'"
            << endl;
    }

    return true;
}


// ************************************************************************* //
