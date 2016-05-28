/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "includeFuncEntry.H"
#include "functionObjectList.H"
#include "dictionary.H"
#include "IFstream.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionEntries::includeFuncEntry::typeName
(
    Foam::functionEntries::includeFuncEntry::typeName_()
);

// Don't lookup the debug switch here as the debug switch dictionary
// might include includeFuncEntry
int Foam::functionEntries::includeFuncEntry::debug(0);

bool Foam::functionEntries::includeFuncEntry::report(false);


namespace Foam
{
namespace functionEntries
{
    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeFuncEntry,
        execute,
        dictionaryIstream
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeFuncEntry,
        execute,
        primitiveEntryIstream
    );
}
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeFuncEntry::funcPath
(
    const word& fName,
    const dictionary& dict
)
{
    // Search the system and etc directories for the file and return the path
    return functionObjectList::findDict(fName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeFuncEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    const word fName(is);
    const fileName fPath(funcPath(fName, parentDict));
    IFstream ifs(fPath);

    if (ifs)
    {
        if (Foam::functionEntries::includeFuncEntry::report)
        {
            Info<< fPath << endl;
        }
        parentDict.read(ifs);
        return true;
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Cannot open functionObject file "
            << (ifs.name().size() ? ifs.name() : fileName(fName))
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}


bool Foam::functionEntries::includeFuncEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    const word fName(is);
    const fileName fPath(funcPath(fName, parentDict));
    IFstream ifs(fPath);

    if (ifs)
    {
        if (Foam::functionEntries::includeFuncEntry::report)
        {
            Info<< fPath << endl;
        }
        entry.read(parentDict, ifs);
        return true;
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Cannot open functionObject file "
            << (ifs.name().size() ? ifs.name() : fileName(fName))
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}


// ************************************************************************* //
