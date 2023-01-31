/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(includeFuncEntry, 0);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeFuncEntry,
        execute,
        dictionaryIstream
    );
}
}


Foam::fileName Foam::functionEntries::includeFuncEntry::functionObjectDictPath
(
    "caseDicts/postProcessing"
);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeFuncEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    // Read line containing the function name and the optional arguments
    const string fNameArgs(readFuncNameArgs(is));

    return readConfigFile
    (
        "function",
        fNameArgs,
        parentDict,
        functionObjectDictPath,
        {"file", is.name() + " at line " + Foam::name(is.lineNumber())}
    );
}


// ************************************************************************* //
