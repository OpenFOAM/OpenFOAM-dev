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

#include "primitiveEntry.H"
#include "dictionaryEntry.H"
#include "functionEntry.H"
#include "includeEntry.H"
#include "inputModeEntry.H"
#include "stringOps.H"

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::entry::getKeyword(keyType& keyword, Istream& is)
{
    token keywordToken;

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

    // If the token is a valid keyword set 'keyword' return true...
    if (keywordToken.isWord())
    {
        keyword = keywordToken.wordToken();
        return true;
    }
    else if (keywordToken.isString())
    {
        // Enable wildcards
        keyword = keywordToken.stringToken();
        return true;
    }
    // If it is the end of the dictionary or file return false...
    else if (keywordToken == token::END_BLOCK || is.eof())
    {
        return false;
    }
    // Otherwise the token is invalid
    else
    {
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


bool Foam::entry::New(dictionary& parentDict, Istream& is)
{
    is.fatalCheck("entry::New(const dictionary& parentDict, Istream&)");

    keyType keyword;

    // Get the next keyword and if invalid return false
    if (!getKeyword(keyword, is))
    {
        return false;
    }
    else  // Keyword starts entry ...
    {
        if
        (
           !disableFunctionEntries
         && keyword[0] == '#'
        )                           // ... Function entry
        {
            word functionName = keyword(1, keyword.size()-1);
            return functionEntry::execute(functionName, parentDict, is);
        }
        else if
        (
           !disableFunctionEntries
         && keyword[0] == '$'
        )                           // ... Substitution entry
        {
            token nextToken(is);
            is.putBack(nextToken);

            if (keyword.size() > 2 && keyword[1] == token::BEGIN_BLOCK)
            {
                // Recursive substitution mode. Replace between {} with
                // expansion and then let standard variable expansion deal
                // with rest.
                string s(keyword(2, keyword.size()-3));
                // Substitute dictionary and environment variables. Do not allow
                // empty substitutions.
                stringOps::inplaceExpand(s, parentDict, true, false);
                keyword.std::string::replace(1, keyword.size()-1, s);
            }

            if (nextToken == token::BEGIN_BLOCK)
            {
                word varName = keyword(1, keyword.size()-1);

                // lookup the variable name in the given dictionary
                const entry* ePtr = parentDict.lookupScopedEntryPtr
                (
                    varName,
                    true,
                    true
                );

                if (ePtr)
                {
                    // Read as primitiveEntry
                    const keyType newKeyword(ePtr->stream());

                    return parentDict.add
                    (
                        new dictionaryEntry(newKeyword, parentDict, is),
                        false
                    );
                }
                else
                {
                    FatalIOErrorIn
                    (
                        "entry::New(const dictionary& parentDict, Istream&)",
                        is
                    )
                        << "Attempt to use undefined variable " << varName
                        << " as keyword"
                        << exit(FatalIOError);
                    return false;
                }
            }
            else
            {
                parentDict.substituteScopedKeyword(keyword);
            }

            return true;
        }
        else if
        (
           !disableFunctionEntries
         && keyword == "include"
        )                           // ... For backward compatibility
        {
            return functionEntries::includeEntry::execute(parentDict, is);
        }
        else                        // ... Data entries
        {
            token nextToken(is);
            is.putBack(nextToken);

            // Deal with duplicate entries
            bool mergeEntry = false;

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
                    // clear dictionary so merge acts like overwrite
                    if (existingPtr->isDict())
                    {
                        existingPtr->dict().clear();
                    }
                    mergeEntry = true;
                }
                else if (functionEntries::inputModeEntry::protect())
                {
                    // read and discard the entry
                    if (nextToken == token::BEGIN_BLOCK)
                    {
                        dictionaryEntry dummy(keyword, parentDict, is);
                    }
                    else
                    {
                        primitiveEntry  dummy(keyword, parentDict, is);
                    }
                    return true;
                }
                else if (functionEntries::inputModeEntry::error())
                {
                    FatalIOErrorIn
                    (
                        "entry::New(const dictionary& parentDict, Istream&)",
                        is
                    )
                        << "ERROR! duplicate entry: " << keyword
                        << exit(FatalIOError);

                    return false;
                }
            }

            if (nextToken == token::BEGIN_BLOCK)
            {
                return parentDict.add
                (
                    new dictionaryEntry(keyword, parentDict, is),
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

    // Get the next keyword and if invalid return false
    if (!getKeyword(keyword, is))
    {
        return autoPtr<entry>(NULL);
    }
    else // Keyword starts entry ...
    {
        token nextToken(is);
        is.putBack(nextToken);

        if (nextToken == token::BEGIN_BLOCK)
        {
            return autoPtr<entry>
            (
                new dictionaryEntry(keyword, dictionary::null, is)
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
