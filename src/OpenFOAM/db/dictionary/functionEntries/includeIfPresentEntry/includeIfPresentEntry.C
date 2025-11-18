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

#include "includeIfPresentEntry.H"
#include "fileOperation.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(includeIfPresentEntry, 0);

    addToRunTimeSelectionTable
    (
        functionEntry,
        includeIfPresentEntry,
        dictionary
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeIfPresentEntry,
        execute,
        primitiveEntryIstream
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::includeIfPresentEntry::includeIfPresentEntry
(
    const dictionary& parentDict,
    Istream& is
)
:
    includeEntry(typeName, parentDict, is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeIfPresentEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    const fileName fName
    (
        includeFileName(is, includeEntry::fName(), contextDict)
    );

    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEntry::log)
        {
            Info<< fName << endl;
        }
        contextDict.read(ifs);
    }

    return true;
}


bool Foam::functionEntries::includeIfPresentEntry::execute
(
    const dictionary& contextDict,
    primitiveEntry& contextEntry,
    Istream& is
)
{
    const fileName fName(includeFileName(is, fileName(is), contextDict));

    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEntry::log)
        {
            Info<< fName << endl;
        }
        contextEntry.read(contextDict, ifs);
    }

    return true;
}

// ************************************************************************* //
