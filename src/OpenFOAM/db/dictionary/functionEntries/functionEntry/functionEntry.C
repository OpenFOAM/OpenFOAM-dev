/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "functionEntry.H"
#include "IOstreams.H"
#include "ISstream.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineMemberFunctionSelectionTable
    (
        functionEntry,
        execute,
        dictionaryIstream
    );

    defineMemberFunctionSelectionTable
    (
        functionEntry,
        execute,
        primitiveEntryIstream
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::token Foam::functionEntry::readLine(Istream& is)
{
    if (isA<Pstream>(is))
    {
        return token(is);
    }
    else
    {
        return token(word(readFuncNameArgs(is)), is.lineNumber());
    }
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::string Foam::functionEntry::readFuncNameArgs(Istream& is)
{
    string fNameArgs;

    // Read the function name with arguments if on the same line
    const token fName(is);

    if (fName.isWord())
    {
        const word& fNameWord = fName.wordToken();

        if (fNameWord.find(token::BEGIN_LIST) != string::npos)
        {
            // If the function name includes a '(' push it back onto the stream
            // and re-read as a list

            ISstream& iss = dynamic_cast<ISstream&>(is);

            for
            (
                string::const_reverse_iterator rit = fNameWord.rbegin();
                rit != fNameWord.rend();
                ++rit
            )
            {
                iss.putback(*rit);
            }

            iss.readList(fNameArgs);
        }
        else
        {
            // Read the next token to check for '('
            // in case the optional arguments start on the next line

            const token nextToken(is);

            if
            (
                nextToken.isPunctuation()
             && nextToken.pToken() == token::BEGIN_LIST
            )
            {
                ISstream& iss = dynamic_cast<ISstream&>(is);

                iss.putback(token::BEGIN_LIST);
                iss.readList(fNameArgs);
                fNameArgs = fNameWord + fNameArgs;
            }
            else
            {
                is.putBack(nextToken);
                fNameArgs = fNameWord;
            }
        }
    }
    else if (fName.isString())
    {
        // If the function name is a string delimit with '"'s
        fNameArgs = '"' + fName.stringToken() + '"';
    }
    else
    {
        // For any other kind of string return for error reporting
        fNameArgs = fName.anyStringToken();
    }

    return fNameArgs;
}


bool Foam::functionEntry::insert
(
    dictionary& parentDict,
    const string& str
)
{
    parentDict.read(IStringStream(str)());
    return true;
}


bool Foam::functionEntry::insert
(
    const dictionary& parentDict,
    primitiveEntry& thisEntry,
    const string& str
)
{
    thisEntry.read(parentDict, IStringStream(str)());
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntry::functionEntry
(
    const word& key,
    const dictionary& dict,
    Istream& is
)
:
    primitiveEntry(key, readLine(is))
{}


// * * * * * * * * * * * * Member Function Selectors * * * * * * * * * * * * //

bool Foam::functionEntry::execute
(
    const word& functionName,
    dictionary& parentDict,
    Istream& is
)
{
    is.fatalCheck
    (
        "functionEntry::execute"
        "(const word& functionName, dictionary& parentDict, Istream&)"
    );

    if (!executedictionaryIstreamMemberFunctionTablePtr_)
    {
        cerr<< "functionEntry::execute"
            << "(const word&, dictionary&, Istream&)"
            << " not yet initialised, function = "
            << functionName.c_str() << std::endl;

        // Return true to keep reading
        return true;
    }

    executedictionaryIstreamMemberFunctionTable::iterator mfIter =
        executedictionaryIstreamMemberFunctionTablePtr_->find(functionName);

    if (mfIter == executedictionaryIstreamMemberFunctionTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown functionEntry '" << functionName
            << "' in " << is.name() << " near line " << is.lineNumber()
            << nl << nl
            << "Valid functionEntries are :" << endl
            << executedictionaryIstreamMemberFunctionTablePtr_->toc()
            << exit(FatalError);
    }

    return mfIter()(parentDict, is);
}


bool Foam::functionEntry::execute
(
    const word& functionName,
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    is.fatalCheck
    (
        "functionEntry::execute"
        "(const word&, const dictionary&, primitiveEntry&, Istream&)"
    );

    if (!executeprimitiveEntryIstreamMemberFunctionTablePtr_)
    {
        cerr<< "functionEntry::execute"
            << "(const word&, const dictionary&, primitiveEntry&, Istream&)"
            << " not yet initialised, function = "
            << functionName.c_str() << std::endl;

        // return true to keep reading anyhow
        return true;
    }

    executeprimitiveEntryIstreamMemberFunctionTable::iterator mfIter =
        executeprimitiveEntryIstreamMemberFunctionTablePtr_->find(functionName);

    if (mfIter == executeprimitiveEntryIstreamMemberFunctionTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown functionEntry '" << functionName
            << "' in " << is.name() << " near line " << is.lineNumber()
            << nl << nl
            << "Valid functionEntries are :" << endl
            << executeprimitiveEntryIstreamMemberFunctionTablePtr_->toc()
            << exit(FatalError);
    }

    return mfIter()(parentDict, entry, is);
}


void Foam::functionEntry::write(Ostream& os) const
{
    writeKeyword(os, keyword());

    if (size() == 1)
    {
        os << operator[](0) << endl;
    }
    else
    {
        FatalIOErrorInFunction(os)
            << "Incorrect number of tokens in functionEntry, "
               "should be a single word."
            << exit(FatalIOError);
    }
}


// ************************************************************************* //
