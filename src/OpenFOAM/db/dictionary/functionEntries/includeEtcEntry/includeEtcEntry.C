/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "includeEtcEntry.H"
#include "etcFiles.H"
#include "stringOps.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(includeEtcEntry, 0);

    addToRunTimeSelectionTable
    (
        functionEntry,
        includeEtcEntry,
        dictionary
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEtcEntry,
        execute,
        primitiveEntryIstream
    );
}
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeEtcEntry::includeFileName
(
    const fileName& dir,
    const fileName& f,
    const dictionary& dict
) const
{
    fileName fName(f);

    // Substitute dictionary and environment variables.
    // Allow empty substitutions.
    stringOps::inplaceExpandEntry(fName, dict, true, true);

    if (fName.empty() || fName.isAbsolute())
    {
        return fName;
    }
    else
    {
        // Search the etc directories for the file
        return findEtcFile(fName);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::includeEtcEntry::includeEtcEntry
(
    const label lineNumber,
    const dictionary& parentDict,
    Istream& is
)
:
    includeEntry(typeName, lineNumber, parentDict, is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEtcEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    return includeEntry::execute(contextDict, is);
}


bool Foam::functionEntries::includeEtcEntry::execute
(
    const dictionary& contextDict,
    primitiveEntry& contextEntry,
    Istream& is
)
{
    return
        includeEtcEntry(is.lineNumber(), contextDict, is)
       .virtualExecute(contextDict, contextEntry, is);
}


// ************************************************************************* //
