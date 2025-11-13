/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
#include "ISstream.H"

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

void Foam::functionEntry::readRestOfArgs(string& fNameArgs, Istream& is)
{
    int listDepth = 1
        + fNameArgs.count(token::BEGIN_LIST) - fNameArgs.count(token::END_LIST);

    if (listDepth != 0)
    {
        ISstream& iss = dynamic_cast<ISstream&>(is);

        char c;

        while (iss.get(c))
        {
            if (c == token::BEGIN_LIST)
            {
                listDepth++;
            }
            else if (c == token::END_LIST)
            {
                listDepth--;

                if (listDepth <= 0)
                {
                    break;
                }
            }

            fNameArgs += c;
        }
    }
    else
    {
        // Remove trailing token::END_LIST
        fNameArgs = fNameArgs(fNameArgs.size() - 1);
    }
}


Foam::tokenList Foam::functionEntry::readFuncNameArgList(Istream& is)
{
    tokenList funcNameArgList;
    const label lineNumber0 = is.lineNumber();
    string::size_type argsStart = string::npos;

    // Read the function name with arguments if on the same line
    const token fName(is);

    if (fName.isWord() || fName.isString())
    {
        if (fName.isString())
        {
            funcNameArgList.append(fName);
        }
        else
        {
            word fNameArgs = fName.wordToken();

            argsStart = fNameArgs.find(token::BEGIN_LIST);

            // If the function name includes a '(' read the rest of the list
            if (argsStart != string::npos)
            {
                const word fName(fNameArgs(0, argsStart));
                funcNameArgList.append(token(fName, lineNumber0));
                funcNameArgList.append(token(token::BEGIN_LIST, lineNumber0));

                string fArgs
                (
                    fNameArgs(argsStart + 1, fNameArgs.size() - argsStart)
                );
                readRestOfArgs(fArgs, is);

                funcNameArgList.append(token(fArgs, lineNumber0));
                funcNameArgList.append(token(token::END_LIST, lineNumber0));
            }
            else
            {
                funcNameArgList.append(token(fNameArgs, lineNumber0));
            }
        }

        if (argsStart == string::npos)
        {
            // Read the next token to check for '('
            // in case the optional arguments start on the next line

            const label fNameLineNumber = is.lineNumber();
            const token nextToken(is);

            if
            (
                nextToken.isPunctuation()
             && nextToken.pToken() == token::BEGIN_LIST
            )
            {
                funcNameArgList.append(token(token::BEGIN_LIST, lineNumber0));

                const token nextToken(is);

                // Check if the argument list has been parsed previously
                // and converted into a single string token
                if (nextToken.isString())
                {
                    funcNameArgList.append(nextToken);

                    const token endListToken(is);
                    if
                    (
                        endListToken.isPunctuation()
                     && endListToken.pToken() == token::END_LIST
                    )
                    {
                        funcNameArgList.append(endListToken);
                    }
                    else
                    {
                        FatalIOErrorInFunction(is)
                            << "Unclosed argument list " << nextToken
                            << " in functionEntry " << funcNameArgList
                            << exit(FatalIOError);
                    }
                }
                else
                {
                    string fArgs(nextToken.anyStringToken());
                    readRestOfArgs(fArgs, is);

                    // Reinstate the \n between fName and the '('
                    if (fNameLineNumber != lineNumber0)
                    {
                        fArgs =
                            string(fNameLineNumber - lineNumber0, '\n')
                          + fArgs;
                    }

                    funcNameArgList.append(token(fArgs, fNameLineNumber));
                    funcNameArgList.append
                    (
                        token(token::END_LIST, is.lineNumber())
                    );
                }
            }
            else
            {
                is.putBack(nextToken);
            }
        }
    }
    else
    {
        // For any other kind of string return for error reporting
        funcNameArgList.append(fName);
    }

    return funcNameArgList;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

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
    const keyType& key,
    const dictionary& dict,
    Istream& is
)
:
    primitiveEntry(key, readFuncNameArgList(is))
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
            << "' in " << is.name() << " at line " << is.lineNumber()
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

        // Return true to keep reading
        return true;
    }

    executeprimitiveEntryIstreamMemberFunctionTable::iterator mfIter =
        executeprimitiveEntryIstreamMemberFunctionTablePtr_->find(functionName);

    if (mfIter == executeprimitiveEntryIstreamMemberFunctionTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown functionEntry '" << functionName
            << "' in " << is.name() << " at line " << is.lineNumber()
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

    if (size())
    {
        primitiveEntry::write(os, true);
    }
    else
    {
        FatalIOErrorInFunction(os)
            << "Empty tokenList in functionEntry"
            << exit(FatalIOError);
    }
}


// ************************************************************************* //
