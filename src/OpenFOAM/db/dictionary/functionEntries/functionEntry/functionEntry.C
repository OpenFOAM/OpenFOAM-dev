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
#include "dummyEntry.H"
#include "ifEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineFunctionTypeName(functionEntry);
    defineRunTimeSelectionTable(functionEntry, dictionary);

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

    if (dynamic_cast<ISstream*>(&is))
    {
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

                    if (listDepth == 0)
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
            fNameArgs.resize(fNameArgs.size() - 1);
        }
    }
    else
    {
        OStringStream argStream;

        token currToken;

        while
        (
            !is.read(currToken).bad()
            && currToken.good()
            && !(currToken == token::END_STATEMENT && listDepth == 1)
        )
        {
            if (currToken == token::BEGIN_LIST)
            {
                listDepth++;
            }
            else if (currToken == token::END_LIST)
            {
                listDepth--;

                if (listDepth == 0)
                {
                    break;
                }
            }

            argStream << currToken;
        }

        fNameArgs.append(argStream.str());
    }
}


Foam::tokenList Foam::functionEntry::readArgList
(
    const functionName& functionType,
    Istream& is,
    const bool optional
)
{
    tokenList argList;

    // Read the next token to check for '('
    // in case the optional arguments start on the next line
    const token nextToken(is);

    if
    (
        nextToken.isPunctuation()
     && nextToken.pToken() == token::BEGIN_LIST
    )
    {
        argList.append(nextToken);

        const token nextToken(is);

        // Check if the argument list has been parsed previously
        // and converted into a single string token
        if (nextToken.isString())
        {
            argList.append(nextToken);

            const token endListToken(is);
            if
            (
                endListToken.isPunctuation()
             && endListToken.pToken() == token::END_LIST
            )
            {
                argList.append(endListToken);
            }
            else
            {
                FatalIOErrorInFunction(is)
                    << "Unclosed argument list " << nextToken
                    << " in functionEntry " << argList
                    << exit(FatalIOError);
            }
        }
        else
        {
            string fArgs(nextToken.anyStringToken());
            readRestOfArgs(fArgs, is);
            argList.append(token(fArgs, nextToken.lineNumber()));
            argList.append(token(token::END_LIST, is.lineNumber()));
        }
    }
    else
    {
        if (!optional)
        {
            FatalIOErrorInFunction(is)
                << "Expected " << char(token::BEGIN_LIST)
                << " to open argument list, but found " << nextToken
                << " in functionEntry " << functionType
                << exit(FatalIOError);
        }

        is.putBack(nextToken);
    }

    return argList;
}


Foam::tokenList Foam::functionEntry::readFuncNameArgList
(
    const functionName& functionType,
    Istream& is
)
{
    tokenList funcNameArgList;
    string::size_type argsStart = string::npos;

    // Read the function name with arguments if on the same line
    const token fName(is);
    const label fNameLineNumber = is.lineNumber();

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
                funcNameArgList.append(token(fName, fNameLineNumber));
                funcNameArgList.append
                (
                    token(token::BEGIN_LIST, fNameLineNumber)
                );

                string fArgs
                (
                    fNameArgs(argsStart + 1, fNameArgs.size() - argsStart)
                );
                readRestOfArgs(fArgs, is);

                funcNameArgList.append(token(fArgs, fNameLineNumber));
                funcNameArgList.append(token(token::END_LIST, fNameLineNumber));
            }
            else
            {
                funcNameArgList.append(token(fNameArgs, fNameLineNumber));
            }
        }

        if (argsStart == string::npos)
        {
            funcNameArgList.append(readArgList(functionType, is, true));
        }
    }
    else
    {
        // For any other kind of string return for error reporting
        funcNameArgList.append(fName);
    }

    return funcNameArgList;
}


Foam::tokenList Foam::functionEntry::readFileNameArgList
(
    const functionName& functionType,
    Istream& is
)
{
    tokenList fileNameArgList;

    // Read the file name
    const token fName(is);

    // Check the file name is a string
    if (fName.isString())
    {
        // Append the file name
        fileNameArgList.append(fName);

        // Append the optional argument list
        fileNameArgList.append(readArgList(functionType, is, true));
    }
    else
    {
        // For any other kind of string return for error reporting
        fileNameArgList.append(fName);
    }

    return fileNameArgList;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::functionEntry::insert
(
    const dictionary& contextDict,
    dictionary& context,
    const token& t,
    Istream& is
)
{
    is.putBack(t);
    return entry::New(context, is);
}


bool Foam::functionEntry::insert
(
    const dictionary& contextDict,
    primitiveEntry& context,
    const token& t,
    Istream& is
)
{
    context.append(t, contextDict, is);
    return true;
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
    const keyType& key,
    const dictionary& dict
)
:
    primitiveEntry(key)
{}

Foam::functionEntry::functionEntry
(
    const keyType& key,
    const dictionary& dict,
    Istream& is
)
:
    primitiveEntry(key, readFuncNameArgList(typeName, is))
{}


Foam::functionEntry::functionEntry
(
    const functionName& functionType,
    const dictionary& dict
)
:
    primitiveEntry(keyType(functionType))
{}


Foam::functionEntry::functionEntry
(
    const functionName& functionType,
    const dictionary& dict,
    const token& token
)
:
    primitiveEntry(keyType(functionType), token)
{}


Foam::functionEntry::functionEntry
(
    const functionName& functionType,
    const dictionary& dict,
    const tokenList& tokens
)
:
    primitiveEntry(keyType(functionType), tokens)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::functionEntry> Foam::functionEntry::New
(
    const keyType& functionName,
    const dictionary& parentDict,
    Istream& is
)
{
    // If the functionEntry constructor table has not yet been constructed
    // ignore functionEntries by instantiating a dummy entry for each
    if (!dictionaryConstructorTablePtr_)
    {
        return autoPtr<functionEntry>
        (
            new functionEntries::dummyEntry(functionName, parentDict, is)
        );
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(functionName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown functionEntry "
            << functionName << nl << nl
            << "Valid functions are : " << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << endl; // exit(FatalError);

        return autoPtr<functionEntry>
        (
            new functionEntry(functionName, parentDict, is)
        );
    }

    return autoPtr<functionEntry>(cstrIter()(parentDict, is));
}


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
}


// ************************************************************************* //
