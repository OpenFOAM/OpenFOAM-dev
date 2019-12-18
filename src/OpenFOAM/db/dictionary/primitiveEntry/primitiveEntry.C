/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
#include "dictionary.H"
#include "OSspecific.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveEntry::append(const UList<token>& varTokens)
{
    forAll(varTokens, i)
    {
        newElmt(tokenIndex()++) = varTokens[i];
    }
}


bool Foam::primitiveEntry::expandVariable
(
    const variable& w,
    const dictionary& dict
)
{
    if (w.size() > 2 && w[0] == '$' && w[1] == token::BEGIN_BLOCK)
    {
        // Recursive substitution mode. Replace between {} with expansion.
        string s(w(2, w.size() - 3));

        // Substitute dictionary and environment variables. Do not allow
        // empty substitutions.
        stringOps::inplaceExpand(s, dict, true, false);
        variable newW(w);
        newW.std::string::replace(1, newW.size() - 1, s);

        return expandVariable(newW, dict);
    }
    else
    {
        string varName = w(1, w.size() - 1);

        // Lookup the variable name in the given dictionary....
        // Note: allow wildcards to match? For now disabled since following
        // would expand internalField to wildcard match and not expected
        // internalField:
        //      internalField XXX;
        //      boundaryField { ".*" {YYY;} movingWall {value $internalField;}
        const entry* ePtr = dict.lookupScopedEntryPtr(varName, true, false);

        // ...if defined append its tokens into this
        if (ePtr)
        {
            if (ePtr->isDict())
            {
                append(ePtr->dict().tokens());
            }
            else
            {
                append(ePtr->stream());
            }
        }
        else
        {
            // Not in the dictionary - try an environment variable
            string envStr = getEnv(varName);

            if (envStr.empty())
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Illegal dictionary entry or environment variable name "
                    << varName << endl << "Valid dictionary entries are "
                    << dict.toc() << exit(FatalIOError);

                return false;
            }
            append(tokenList(IStringStream('(' + envStr + ')')()));
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::primitiveEntry::primitiveEntry(const keyType& key, const ITstream& is)
:
    entry(key),
    ITstream(is)
{
    name() += '/' + keyword();
}


Foam::primitiveEntry::primitiveEntry(const keyType& key, const token& t)
:
    entry(key),
    ITstream(key, tokenList(1, t))
{}


Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    const UList<token>& tokens
)
:
    entry(key),
    ITstream(key, tokens)
{}


Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    List<token>&& tokens
)
:
    entry(key),
    ITstream(key, move(tokens))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::primitiveEntry::startLineNumber() const
{
    const tokenList& tokens = *this;

    if (tokens.empty())
    {
        return -1;
    }
    else
    {
        return tokens.first().lineNumber();
    }
}


Foam::label Foam::primitiveEntry::endLineNumber() const
{
    const tokenList& tokens = *this;

    if (tokens.empty())
    {
        return -1;
    }
    else
    {
        return tokens.last().lineNumber();
    }
}


Foam::ITstream& Foam::primitiveEntry::stream() const
{
    ITstream& is = const_cast<primitiveEntry&>(*this);
    is.rewind();
    return is;
}


const Foam::dictionary& Foam::primitiveEntry::dict() const
{
    FatalErrorInFunction
        << "Attempt to return primitive entry " << info()
        << " as a sub-dictionary"
        << abort(FatalError);

    return dictionary::null;
}


Foam::dictionary& Foam::primitiveEntry::dict()
{
    FatalErrorInFunction
        << "Attempt to return primitive entry " << info()
        << " as a sub-dictionary"
        << abort(FatalError);

    return const_cast<dictionary&>(dictionary::null);
}


// ************************************************************************* //
