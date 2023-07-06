/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
#include "calcIncludeEntry.H"
#include "dictionary.H"
#include "dynamicCode.H"
#include "codeStream.H"
#include "addToMemberFunctionSelectionTable.H"

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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::string Foam::functionEntries::calcEntry::calc
(
    const dictionary& dict,
    Istream& is
)
{
    Info<< "Expanding #calc at line " << is.lineNumber()
        << " in file " <<  dict.name() << endl;

    dynamicCode::checkSecurity
    (
        "functionEntries::calcEntry::execute(..)",
        dict
    );

    // Construct codeDict for codeStream
    // with dict as parent dictionary for string expansion
    // and variable substitution
    dictionary codeDict(dict, dictionary());

    // Read the code expression string delimited by either '"..."' or '#{...#}'
    token t(is);

    if (t.isVerbatimString())
    {
        const verbatimString& s = t.verbatimStringToken();

        calcIncludeEntry::codeInclude(codeDict);
        codeDict.add("code", s);
    }
    else if (t.isString())
    {
        const string& s = t.stringToken();

        calcIncludeEntry::codeInclude(codeDict);
        codeDict.add("code", "os << (" + s + ");");
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong string type for #calc" << nl
            << "    Expected either a string delimited by '\"...\"' "
               "or a verbatim string delimited by '#{...#}' " << nl
            << "    found token " << t
            << exit(FatalIOError);
    }

    codeStream::streamingFunctionType function = codeStream::getFunction
    (
        dict,
        codeDict
    );

    // Use function to write stream
    OStringStream os(is.format());
    (*function)(os, dict);

    // Return the string containing the results of the calculation
    return os.str();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::calcEntry::execute
(
    dictionary& dict,
    Istream& is
)
{
    return insert(dict, calc(dict, is));
}


bool Foam::functionEntries::calcEntry::execute
(
    const dictionary& dict,
    primitiveEntry& thisEntry,
    Istream& is
)
{
    return insert(dict, thisEntry, calc(dict, is));
}


// ************************************************************************* //
