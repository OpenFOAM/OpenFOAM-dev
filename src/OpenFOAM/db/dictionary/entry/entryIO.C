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

#include "entry.H"
#include "dictionaryListEntry.H"
#include "includeEntry.H"
#include "inputModeEntry.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::entry::getKeyword
(
    keyType& keyword,
    token& keywordToken,
    label& keywordLineNo,
    Istream& is
)
{
    // Read the next valid token discarding spurious ';'s
    do
    {
        if
        (
            is.read(keywordToken).bad()
         || is.eof()
         || !keywordToken.good()
        )
        {
            return false;
        }
    }
    while (keywordToken == token::END_STATEMENT);

    keyword = keywordToken;
    keywordLineNo = is.lineNumber();

    return !keyword.isUndefined();
}


bool Foam::entry::getKeyword
(
    keyType& keyword,
    label& keywordLineNo,
    Istream& is
)
{
    token keywordToken;
    bool ok = getKeyword(keyword, keywordToken, keywordLineNo, is);

    if (ok)
    {
        return true;
    }
    else
    {
        // Do some more checking
        if (keywordToken == token::END_BLOCK || is.eof())
        {
            return false;
        }
        else
        {
            // Otherwise the token is invalid
            cerr<< "--> FOAM Warning : " << std::endl
                << "    From function "
                << "entry::getKeyword(keyType&, Istream&)" << std::endl
                << "    in file " << __FILE__
                << " at line " << __LINE__ << std::endl
                << "    Reading " << is.name().c_str() << std::endl
                << "    found " << keywordToken << std::endl
                << "    expected either " << token::END_BLOCK << " or EOF"
                << std::endl;
            return false;
        }
    }
}


bool Foam::entry::New(dictionary& parentDict, Istream& is)
{
    is.fatalCheck("entry::New(const dictionary& parentDict, Istream&)");

    keyType keyword;
    label keywordLineNo = -1;
    token keyToken;

    // Get the next keyword and if valid return true
    bool valid = getKeyword(keyword, keyToken, keywordLineNo, is);

    if (!valid)
    {
        // Do some more checking
        if (keyToken == token::END_BLOCK || is.eof())
        {
            return false;
        }
        else if
        (
            keyToken.isLabel()
         || (keyToken.isPunctuation() && keyToken.pToken() == token::BEGIN_LIST)
        )
        {
            is.putBack(keyToken);
            return parentDict.add
            (
                new dictionaryListEntry(parentDict, is),
                false
            );
        }
        else
        {
            // Otherwise the token is invalid
            cerr<< "--> FOAM Warning : " << std::endl
                << "    From function "
                << "entry::New(dictionary&, Istream&)" << std::endl
                << "    in file " << __FILE__
                << " at line " << __LINE__ << std::endl
                << "    Reading " << is.name().c_str() << std::endl
                << "    found " << keyToken << std::endl
                << "    expected either " << token::END_BLOCK << " or EOF"
                << std::endl;
            return false;
        }
    }
    else  // Keyword starts entry ...
    {
        if (keyword.isFunctionName())      // ... Function entry
        {
            const word functionName = keyword(1, keyword.size() - 1);

            if (disableFunctionEntries)
            {
                bool success = parentDict.add
                (
                    new functionEntry
                    (
                        keyword,
                        parentDict,
                        is
                    ),
                    false
                );

                return success;
            }
            else
            {
                return functionEntry::execute(functionName, parentDict, is);
            }
        }
        else if
        (
           !disableFunctionEntries
         && keyword.isVariable()
        )                           // ... Substitution entry
        {
            token nextToken(is);
            is.putBack(nextToken);

            if (keyword.size() > 2 && keyword[1] == token::BEGIN_BLOCK)
            {
                // Recursive substitution mode. Replace between {} with
                // expansion and then let standard variable expansion deal
                // with rest.
                string s(keyword(2, keyword.size() - 3));

                // Substitute dictionary and environment variables. Do not allow
                // empty substitutions.
                stringOps::inplaceExpandEntry(s, parentDict, true, false);
                keyword.std::string::replace(1, keyword.size() - 1, s);
            }

            if (nextToken == token::BEGIN_BLOCK)
            {
                const word varName = keyword(1, keyword.size() - 1);

                // Lookup the variable name in the given dictionary
                const entry* ePtr = parentDict.lookupScopedEntryPtr
                (
                    varName,
                    true,
                    true
                );

                if (ePtr == nullptr)
                {
                    FatalIOErrorInFunction(is)
                        << "Attempt to use undefined variable " << varName
                        << " as keyword"
                        << exit(FatalIOError);
                    return false;
                }
                else if (ePtr->isStream())
                {
                    // Read as primitiveEntry
                    const keyType newKeyword(ePtr->stream());

                    return parentDict.add
                    (
                        new dictionaryEntry(newKeyword, parentDict, is),
                        false
                    );
                }
                else if (ePtr->isDict())
                {
                    FatalIOErrorInFunction(is)
                        << "Attempt to substitute sub-dictionary"
                        << nl << *ePtr
                        << "for keyword " << varName
                        << exit(FatalIOError);
                    return false;
                }
            }
            else
            {
                const word varName = keyword(1, keyword.size() - 1);

                // Lookup the variable name in the given dictionary
                const entry* ePtr = parentDict.lookupScopedEntryPtr
                (
                    varName,
                    true,
                    true
                );

                if (ePtr == nullptr)
                {
                    FatalIOErrorInFunction(is)
                        << "Attempt to use undefined variable " << keyword
                        << " as a keyword"
                        << exit(FatalIOError);
                    return false;
                }
                else if (ePtr->isDict())
                {
                    // Insert the sub-dict entries into this dictionary
                    const dictionary& addDict = ePtr->dict();

                    forAllConstIter(IDLList<entry>, addDict, iter)
                    {
                        parentDict.add(iter());
                    }
                    return true;
                }
                else
                {
                    FatalIOErrorInFunction(is)
                        << "Attempt to substitute primitive entry "
                        << *ePtr << "as a sub-dictionary"
                        << exit(FatalIOError);
                    return false;
                }
            }

            return true;
        }
        else                        // ... Data entries
        {
            token nextToken(is);
            is.putBack(nextToken);

            // Deal with duplicate entries
            bool mergeEntry = false;

            // If function entries are disabled allow duplicate entries
            if (disableFunctionEntries)
            {
                mergeEntry = false;
            }
            else
            {
                // See (using exact match) if entry already present
                entry* existingPtr = parentDict.lookupEntryPtr
                (
                    keyword,
                    false,
                    false
                );

                if (existingPtr)
                {
                    if (functionEntries::inputModeEntry::merge())
                    {
                        mergeEntry = true;
                    }
                    else if (functionEntries::inputModeEntry::overwrite())
                    {
                        // Clear dictionary so merge acts like overwrite
                        if (existingPtr->isDict())
                        {
                            existingPtr->dict().clear();
                        }
                        mergeEntry = true;
                    }
                    else if (functionEntries::inputModeEntry::protect())
                    {
                        // Read and discard the entry
                        if (nextToken == token::BEGIN_BLOCK)
                        {
                            dictionaryEntry dummy
                            (
                                keyword,
                                keywordLineNo,
                                parentDict,
                                is
                            );
                        }
                        else
                        {
                            primitiveEntry dummy(keyword, parentDict, is);
                        }
                        return true;
                    }
                    else if (functionEntries::inputModeEntry::error())
                    {
                        FatalIOErrorInFunction(is)
                            << "ERROR! duplicate entry: " << keyword
                            << exit(FatalIOError);
                        return false;
                    }
                }
            }

            if (nextToken == token::BEGIN_BLOCK)
            {
                return parentDict.add
                (
                    new dictionaryEntry(keyword, keywordLineNo, parentDict, is),
                    mergeEntry
                );
            }
            else
            {
                return parentDict.add
                (
                    new primitiveEntry(keyword, parentDict, is),
                    mergeEntry
                );
            }
        }
    }
}


Foam::autoPtr<Foam::entry> Foam::entry::New(Istream& is)
{
    is.fatalCheck("entry::New(Istream&)");

    keyType keyword;
    label keywordLineNo = -1;

    // Get the next keyword and if invalid return false
    if (!getKeyword(keyword, keywordLineNo, is))
    {
        return autoPtr<entry>(nullptr);
    }
    else // Keyword starts entry ...
    {
        token nextToken(is);
        is.putBack(nextToken);

        if (nextToken == token::BEGIN_BLOCK)
        {
            return autoPtr<entry>
            (
                new dictionaryEntry
                (
                    keyword,
                    keywordLineNo,
                    dictionary::null,
                    is
                )
            );
        }
        else
        {
            return autoPtr<entry>
            (
                new primitiveEntry(keyword, is)
            );
        }
    }
}


// * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const entry& e)
{
    e.write(os);
    return os;
}


// ************************************************************************* //
