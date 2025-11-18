/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(includeFuncEntry, 0);
    addToRunTimeSelectionTable(functionEntry, includeFuncEntry, dictionary);
}
}


Foam::fileName Foam::functionEntries::includeFuncEntry::functionObjectDictPath
(
    "caseDicts/functions"
);


Foam::fileName
Foam::functionEntries::includeFuncEntry::functionObjectTemplatePath
(
    "caseDicts/functionTemplates"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::includeFuncEntry::includeFuncEntry
(
    const functionName& functionType,
    const dictionary& parentDict,
    Istream& is
)
:
    functionEntry
    (
        functionType,
        parentDict,
        is,
        readFuncNameArgList(functionType, is)
    )
{}


Foam::functionEntries::includeFuncEntry::includeFuncEntry
(
    const dictionary& parentDict,
    Istream& is
)
:
    includeFuncEntry(typeName, parentDict, is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeFuncEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    return readConfigFile
    (
        "function",
        // Read line containing the function name and the optional arguments
        funcNameArgs(),
        contextDict,
        functionObjectDictPath,
        "system"
    );
}


// ************************************************************************* //
