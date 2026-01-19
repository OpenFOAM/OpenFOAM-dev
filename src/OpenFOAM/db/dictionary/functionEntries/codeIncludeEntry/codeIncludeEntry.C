/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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

#include "codeIncludeEntry.H"
#include "stringOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(codeIncludeEntry, 0);
    addToRunTimeSelectionTable(functionEntry, codeIncludeEntry, dictionary);

    addBackwardCompatibleToRunTimeSelectionTable
    (
        functionEntry,
        codeIncludeEntry,
        dictionary,
        calcIncludeEntry,
        "#calcInclude"
    );
}
}

// Construct the static include file name cache
Foam::DynamicList<Foam::fileName>
    Foam::functionEntries::codeIncludeEntry::includeFiles_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionEntries::codeIncludeEntry::appendFileName
(
    const dictionary& contextDict,
    const fileName& fName
) const
{
    // Copy the file name for inplace expansion
    fileName expandedFname(fName);

    // Substitute dictionary and environment variables. Allow empty
    // substitutions.
    stringOps::inplaceExpandEntry(expandedFname, contextDict, true, true);

    // Add the file name to the cache
    includeFiles_.append(expandedFname);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::codeIncludeEntry::codeIncludeEntry
(
    const label lineNumber,
    const dictionary& parentDict,
    Istream& is
)
:
    functionEntry
    (
        typeName,
        lineNumber,
        parentDict,
        is,
        readArgOrList(typeName, is)
    ),
    fNames_(readList<fileName>(stream()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::codeIncludeEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    forAll(fNames_, i)
    {
        appendFileName(contextDict, fNames_[i]);
    }

    return true;
}


void Foam::functionEntries::codeIncludeEntry::clear()
{
    includeFiles_.clear();
}


void Foam::functionEntries::codeIncludeEntry::codeInclude(dictionary& codeDict)
{
    if (includeFiles_.size())
    {
        verbatimString codeInclude;
        forAll(includeFiles_, i)
        {
            codeInclude += "#include \"" + includeFiles_[i] + '"' + '\n';
        }

        codeDict.add("codeInclude", codeInclude);
    }
}


// ************************************************************************* //
