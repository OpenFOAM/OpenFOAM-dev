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

#include "codeBlockEntry.H"
#include "endCodeBlockEntry.H"
#include "codeBlockCalcEntry.H"
#include "codeStream.H"
#include "streamEntry.H"
#include "calcEntry.H"
#include "codeIncludeEntry.H"
#include "negEntry.H"
#include "OTstream.H"
#include "ITstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(codeBlockEntry, 0);
    addToRunTimeSelectionTable(functionEntry, codeBlockEntry, dictionary);
}
}

const Foam::word Foam::functionEntries::codeBlockEntry::codeTemplateC =
    "codeBlockTemplate.C";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::codeBlockEntry::codeBlockEntry
(
    const label lineNumber,
    const dictionary& parentDict,
    Istream& is
)
:
    functionEntry(typeName, lineNumber, parentDict),
    codeBlockName_(word::null),
    lib_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::codeBlockEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    if (debug)
    {
        Info<< "Compiling code within " << typeName << " starting at line "
            << is.lineNumber() << " in file " <<  contextDict.name() << endl;
    }

    // Construct codeDict for codeStream with the parent dictionary provided for
    // string expansion and variable substitution and the same name as the
    // parent for consistent error messaging
    dictionary codeDict(typeName, contextDict);
    string codeString;

    // Disable functionEntry expansion
    entry::disableFunctionEntries = true;

    // Read the dictionary up to #endCodeBlock
    OTstream dictStream(contextDict.name());
    label codeIndex = 0;

    token t;
    while (!is.eof() && !is.read(t).bad() && t.good())
    {
        if (t.isFunctionName())
        {
            if (t.functionNameToken() == endCodeBlockEntry::typeName)
            {
                break;
            }
            else if (t.functionNameToken() == codeStream::typeName)
            {
                // Accumulate the #codeStream code strings
                // into a single code block
                codeString +=
                    codeStream::codeString(codeIndex, contextDict, is);

                // Replace the #codeStream with a #codeBlockCalc followed by the
                // pointer to this codeBlockEntry and the index of #codeStream
                dictStream
                    << keyType(codeBlockCalcEntry::typeName) << token::SPACE
                    << reinterpret_cast<uint64_t>(this) << token::SPACE
                    << codeIndex++ << endl;
            }
            else if (t.functionNameToken() == streamEntry::typeName)
            {
                // Accumulate the #stream code strings into a single code block
                codeString +=
                    streamEntry::codeString(codeIndex, codeDict, is);

                // Replace the #stream with a #codeBlockCalc followed by
                // the pointer to this codeBlockEntry and the index of #stream
                dictStream
                    << keyType(codeBlockCalcEntry::typeName) << token::SPACE
                    << reinterpret_cast<uint64_t>(this) << token::SPACE
                    << codeIndex++ << endl;
            }
            else if (t.functionNameToken() == calcEntry::typeName)
            {
                // Accumulate the #calc code strings into a single code block
                codeString +=
                    calcEntry::codeString(codeIndex, codeDict, is);

                // Replace the #calc with a #codeBlockCalc followed by
                // the pointer to this codeBlockEntry and the index of #calc
                dictStream
                    << keyType(codeBlockCalcEntry::typeName) << token::SPACE
                    << reinterpret_cast<uint64_t>(this) << token::SPACE
                    << codeIndex++ << endl;
            }
            else if (t.functionNameToken() == codeIncludeEntry::typeName)
            {
                codeIncludeEntry(t.lineNumber(), contextDict, is).execute
                (
                    contextDict,
                    is
                );
            }
            else if (t.functionNameToken() == negEntry::typeName)
            {
                dictStream.append(t);
            }
            else
            {
                FatalIOErrorInFunction(is)
                    << t.functionNameToken() << " not supported within "
                    << codeBlockEntry::typeName << "..."
                    << endCodeBlockEntry::typeName
                    << exit(FatalIOError);
            }
        }
        else
        {
            dictStream.append(t);
        }
    }

    // Read the dictionary without expanding functionEntries
    // to provide entries for the code to lookup
    contextDict.read(ITstream(dictStream.name(), dictStream)());

    // Add any optional include statements to codeDict
    codeIncludeEntry::codeInclude(codeDict);

    // Add the code entry to the codeDict
    codeDict.add(primitiveEntry("code", codeString, 0));

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

    // Compile the codeBlock library and cache the pointer
    lib_ = codeStream::compile
    (
        contextDict,
        codeDict,
        codeTemplateC,
        codeBlockName_
    );

    // Re-enable functionEntry expansion
    entry::disableFunctionEntries = false;

    // Re-read the dictionary expanding all functionEntries
    // including #codeBlockCalc and #codeBlockStream using the functions in lib_
    contextDict.read(ITstream(dictStream.name(), dictStream)());

    return true;
}


// ************************************************************************* //
