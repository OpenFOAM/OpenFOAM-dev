/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "calcIncludeEntry.H"
#include "dictionary.H"
#include "IFstream.H"
#include "addToMemberFunctionSelectionTable.H"
#include "stringOps.H"
#include "IOobject.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(calcIncludeEntry, 0);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        calcIncludeEntry,
        execute,
        dictionaryIstream
    );
}
}

// Construct the static include file name cache
Foam::DynamicList<Foam::fileName>
    Foam::functionEntries::calcIncludeEntry::includeFiles_;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::calcIncludeEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    // Read the include file name
    fileName fName(is);

    // Substitute dictionary and environment variables. Allow empty
    // substitutions.
    stringOps::inplaceExpandEntry(fName, parentDict, true, true);

    // Add the file name to the cache
    includeFiles_.append(fName);

    return true;
}


void Foam::functionEntries::calcIncludeEntry::clear()
{
    includeFiles_.clear();
}


void Foam::functionEntries::calcIncludeEntry::codeInclude(dictionary& codeDict)
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
