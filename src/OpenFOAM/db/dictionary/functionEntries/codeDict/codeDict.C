/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "codeDict.H"
#include "dynamicCode.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(codeDict, 0);

    addToRunTimeSelectionTable(functionEntry, codeDict, dictionary);
}
}

const Foam::word Foam::functionEntries::codeDict::codeTemplateC =
    "codeDictTemplate.C";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::string Foam::functionEntries::codeDict::codeString
(
    const label index,
    const dictionary& contextDict,
    Istream& is
)
{
    return codeStream::codeString
    (
        typeName,
        "CODE_BLOCK_DICT_FUNCTION",
        index,
        contextDict,
        is
    );
}


Foam::functionEntries::codeDict::streamingFunctionType
Foam::functionEntries::codeDict::getFunction
(
    const dictionary& contextDict,
    const dictionary& codeDict
)
{
    word codeName;
    void* lib = compile
    (
        typeName,
        contextDict,
        codeDict,
        codeTemplateC,
        codeName
    );

    // Find the function handle in the library
    const streamingFunctionType function =
        reinterpret_cast<streamingFunctionType>
        (
            dlSym(lib, codeName)
        );

    if (!function)
    {
        FatalIOErrorInFunction
        (
            contextDict
        )   << "Failed looking up symbol " << codeName
            << " in library " << lib << exit(FatalIOError);
    }

    return function;
}


bool Foam::functionEntries::codeDict::resultStream
(
    dictionary& contextDict,
    Istream& is
)
{
    if (debug)
    {
        Info<< "Using " << typeName << " at line " << is.lineNumber()
            << " in file " <<  contextDict.name() << endl;
    }

    dynamicCode::checkSecurity
    (
        "functionEntries::codeDict::execute(..)",
        contextDict
    );

    // Construct codeDict for codeDict using the context dictionary
    // for string expansion and variable substitution
    const dictionary codeDict(typeName, contextDict, is);

    // Compile and link the code library and get the function pointer
    const streamingFunctionType function = getFunction(contextDict, codeDict);

    // Use function to append to contextDict
    (*function)(contextDict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::codeDict::codeDict
(
    const label lineNumber,
    const dictionary& parentDict,
    Istream& is
)
:
    codeStream(typeName, lineNumber, parentDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::codeDict::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    return resultStream(contextDict, is);
}


// ************************************************************************* //
