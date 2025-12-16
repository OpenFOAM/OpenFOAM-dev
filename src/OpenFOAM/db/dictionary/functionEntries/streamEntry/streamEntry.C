/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "streamEntry.H"
#include "codeIncludeEntry.H"
#include "dictionary.H"
#include "dynamicCode.H"
#include "codeStream.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(streamEntry, 0);

    addToRunTimeSelectionTable(functionEntry, streamEntry, dictionary);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        streamEntry,
        execute,
        primitiveEntryIstream
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::string Foam::functionEntries::streamEntry::codeString
(
    const label index,
    const dictionary& codeDict,
    Istream& is,
    const string& startString,
    const string& endString
)
{
    // Read the code expression string delimited by either '"..."' or '#{...#}'
    token t(is);

    if (t.isString() || t.isVerbatimString())
    {
        return
        (
            "CODE_BLOCK_FUNCTION(" + Foam::name(index) + ")\n"
            "{\n"
            "    #line " + Foam::name(t.lineNumber())
               + " \"" + codeDict.name() + "\"\n"
               + startString + t.anyStringToken() + endString +
            "\n}\n\n"
        );
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong string type for " << typeName << nl
            << "    Expected either a string delimited by '\"...\"' "
               "or a verbatim string delimited by '#{...#}' " << nl
            << "    found token " << t
            << exit(FatalIOError);
        return string::null;
    }
}


Foam::OTstream Foam::functionEntries::streamEntry::resultStream
(
    const dictionary& dict,
    Istream& is,
    const string& startString,
    const string& endString
)
{
    if (debug)
    {
        Info<< "Expanding " << typeName << " at line " << is.lineNumber()
            << " in file " <<  dict.name() << endl;
    }

    dynamicCode::checkSecurity
    (
        "functionEntries::streamEntry::execute(..)",
        dict
    );

    // Construct codeDict for codeStream with the parent dictionary provided for
    // string expansion and variable substitution and the same name as the
    // parent for consistent error messaging
    dictionary codeDict(fileName::null, dict);

    // Read the code expression string delimited by either '"..."' or '#{...#}'
    token t(is);

    if (t.isString() || t.isVerbatimString())
    {
        codeIncludeEntry::codeInclude(codeDict);
        codeDict.add
        (
            primitiveEntry
            (
                "code",
                startString + t.anyStringToken() + endString,
                t.lineNumber()
            )
        );
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong string type for " << typeName << nl
            << "    Expected either a string delimited by '\"...\"' "
               "or a verbatim string delimited by '#{...#}' " << nl
            << "    found token " << t
            << exit(FatalIOError);
    }

    // Add compilation options to simplify compilation error messages
    codeDict.add
    (
        primitiveEntry
        (
            "codeOptions",
            "#{ -fno-show-column -fno-diagnostics-show-caret #}",
            0
        )
    );

    codeStream::streamingFunctionType function = codeStream::getFunction
    (
        dict,
        codeDict
    );

    // Use function to write stream
    OTstream ots(is.name(), is.format());
    ots.lineNumber() = is.lineNumber();
    (*function)(ots, dict);

    // Return the OTstream containing the results of the calculation
    return ots;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::streamEntry::streamEntry
(
    const label lineNumber,
    const dictionary& parentDict,
    Istream& is
)
:
    functionEntry(typeName, lineNumber, parentDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::streamEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    return insert(contextDict, resultStream(contextDict, is));
}


bool Foam::functionEntries::streamEntry::execute
(
    const dictionary& contextDict,
    primitiveEntry& contextEntry,
    Istream& is
)
{
    return insert(contextDict, contextEntry, resultStream(contextDict, is));
}


// ************************************************************************* //
