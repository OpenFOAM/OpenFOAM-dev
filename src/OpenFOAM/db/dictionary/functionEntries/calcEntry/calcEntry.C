/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "calcEntry.H"
#include "addToMemberFunctionSelectionTable.H"
#include "dictionary.H"
#include "dynamicCode.H"
#include "codeStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(calcEntry, 0);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        calcEntry,
        execute,
        dictionaryIstream
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        calcEntry,
        execute,
        primitiveEntryIstream
    );

}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::calcEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& thisEntry,
    Istream& is
)
{
    Info<< "Using #calcEntry at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << endl;

    dynamicCode::checkSecurity
    (
        "functionEntries::calcEntry::execute(..)",
        parentDict
    );

    // Read string
    string s(is);
    // Make sure we stop this entry
    //is.putBack(token(token::END_STATEMENT, is.lineNumber()));

    // Construct codeDict for codeStream
    // must reference parent for stringOps::expand to work nicely.
    dictionary codeSubDict;
    codeSubDict.add("code", "os << (" + s + ");");
    dictionary codeDict(parentDict, codeSubDict);

    codeStream::streamingFunctionType function = codeStream::getFunction
    (
        parentDict,
        codeDict
    );

    // use function to write stream
    OStringStream os(is.format());
    (*function)(os, parentDict);

    // get the entry from this stream
    IStringStream resultStream(os.str());
    thisEntry.read(parentDict, resultStream);

    return true;
}


bool Foam::functionEntries::calcEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    Info<< "Using #calcEntry at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << endl;

    dynamicCode::checkSecurity
    (
        "functionEntries::calcEntry::execute(..)",
        parentDict
    );

    // Read string
    string s(is);
    // Make sure we stop this entry
    //is.putBack(token(token::END_STATEMENT, is.lineNumber()));

    // Construct codeDict for codeStream
    // must reference parent for stringOps::expand to work nicely.
    dictionary codeSubDict;
    codeSubDict.add("code", "os << (" + s + ");");
    dictionary codeDict(parentDict, codeSubDict);

    codeStream::streamingFunctionType function = codeStream::getFunction
    (
        parentDict,
        codeDict
    );

    // use function to write stream
    OStringStream os(is.format());
    (*function)(os, parentDict);

    // get the entry from this stream
    IStringStream resultStream(os.str());
    parentDict.read(resultStream);

    return true;
}


// ************************************************************************* //
