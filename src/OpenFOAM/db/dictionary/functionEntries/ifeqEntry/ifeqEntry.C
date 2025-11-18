/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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

#include "ifeqEntry.H"
#include "ifEntry.H"
#include "elifEntry.H"
#include "elseEntry.H"
#include "endifEntry.H"
#include "stringOps.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(ifeqEntry, 0);
    addToRunTimeSelectionTable(functionEntry, ifeqEntry, dictionary);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        ifeqEntry,
        execute,
        primitiveEntryIstream
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tokenList Foam::functionEntries::ifeqEntry::readArgList
(
    Istream& is
) const
{
    tokenList args;

    {
        const token startArgs(is);

        if
        (
            startArgs.isPunctuation()
         && startArgs.pToken() == token::BEGIN_LIST
        )
        {
            args.append(startArgs);
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Expected " << token::BEGIN_LIST
                << " to start the argument list but found " << startArgs
                << " in " << typeName
                << exit(FatalIOError);
        }
    }

    // Read the two arguments
    args.append(token(is));
    args.append(token(is));

    {
        const token endArgs(is);

        if
        (
            endArgs.isPunctuation()
        && endArgs.pToken() == token::END_LIST
        )
        {
            args.append(endArgs);
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Expected " << token::END_LIST
                << " to end the argument list but found " << endArgs
                << " in " << typeName
                << exit(FatalIOError);
        }
    }

    return args;
}


void Foam::functionEntries::ifeqEntry::readToken(token& t, Istream& is) const
{
    // Skip dummy tokens - avoids entry::getKeyword consuming #else, #endif
    do
    {
        if( is.read(t).bad() || is.eof() || !t.good())
        {
            return;
        }
    }
    while (t == token::END_STATEMENT);
}


Foam::token Foam::functionEntries::ifeqEntry::expand
(
    const dictionary& dict,
    const token& t
) const
{
    if (t.isVariable())
    {
        const variable& var = t.variableToken();

        word varName = var(1, var.size() - 1);

        // Lookup the variable name in the given dictionary
        const entry* ePtr = dict.lookupScopedEntryPtr
        (
            varName,
            true,
            true
        );

        if (ePtr)
        {
            return token(ePtr->stream());
        }
        else
        {
            // String expansion. Allow unset variables
            string expanded(var);
            stringOps::inplaceExpandEntry(expanded, dict, true, true);

            // Re-form as a string token so we can compare to string
            return token(expanded, t.lineNumber());
        }
    }
    else
    {
        return t;
    }
}


bool Foam::functionEntries::ifeqEntry::equalToken
(
    const token& t1,
    const token& t2
) const
{
    const bool eqType = (t1.type() == t2.type());

    switch (t1.type())
    {
        case token::UNDEFINED:
            return eqType;

        case token::PUNCTUATION:
            return (eqType && t1.pToken() == t2.pToken());

        case token::WORD:
        case token::FUNCTIONNAME:
        case token::STRING:
        case token::VERBATIMSTRING:
            if (t2.isAnyString())
            {
                return t1.anyStringToken() == t2.anyStringToken();
            }
            else
            {
                return false;
            }

        case token::VARIABLE:
            FatalErrorInFunction
                << "Attempt to compare an un-expanded variable"
                << InfoProxy<token>(t1)
                << exit(FatalIOError);
            return false;

        case token::INTEGER_32:
            if (eqType)
            {
                return t1.integer32Token() == t2.integer32Token();
            }
            else if (t2.isLabel())
            {
                return t1.labelToken() == t2.labelToken();
            }
            else if (t2.isScalar())
            {
                return t1.labelToken() == t2.scalarToken();
            }
            else
            {
                return false;
            }

        case token::INTEGER_64:
            if (eqType)
            {
                return t1.integer64Token() == t2.integer64Token();
            }
            else if (t2.isLabel())
            {
                return t1.labelToken() == t2.labelToken();
            }
            else if (t2.isScalar())
            {
                return t1.labelToken() == t2.scalarToken();
            }
            else
            {
                return false;
            }

        case token::UNSIGNED_INTEGER_32:
            if (eqType)
            {
                return
                    t1.unsignedInteger32Token() == t2.unsignedInteger32Token();
            }
            else if (t2.isLabel())
            {
                return t1.labelToken() == t2.labelToken();
            }
            else if (t2.isScalar())
            {
                return t1.labelToken() == t2.scalarToken();
            }
            else
            {
                return false;
            }

        case token::UNSIGNED_INTEGER_64:
            if (eqType)
            {
                return
                    t1.unsignedInteger64Token() == t2.unsignedInteger64Token();
            }
            else if (t2.isLabel())
            {
                return t1.labelToken() == t2.labelToken();
            }
            else if (t2.isScalar())
            {
                return t1.labelToken() == t2.scalarToken();
            }
            else
            {
                return false;
            }

        case token::FLOAT_SCALAR:
            if (eqType)
            {
                return equal(t1.floatScalarToken(), t2.floatScalarToken());
            }
            else if (t2.isScalar())
            {
                return t1.scalarToken() == t2.scalarToken();
            }
            else
            {
                return false;
            }

        case token::DOUBLE_SCALAR:
            if (eqType)
            {
                return equal(t1.doubleScalarToken(), t2.doubleScalarToken());
            }
            else if (t2.isScalar())
            {
                return t1.scalarToken() == t2.scalarToken();
            }
            else
            {
                return false;
            }

        case token::LONG_DOUBLE_SCALAR:
            if (eqType)
            {
                return equal
                (
                    t1.longDoubleScalarToken(),
                    t2.longDoubleScalarToken()
                );
            }
            else if (t2.isScalar())
            {
                return t1.scalarToken() == t2.scalarToken();
            }
            else
            {
                return false;
            }

        case token::COMPOUND:
            return false;

        case token::ERROR:
            return eqType;
    }

    return false;
}


void Foam::functionEntries::ifeqEntry::skipUntil
(
    DynamicList<filePos>& stack,
    const dictionary& contextDict,
    const functionName& endWord,
    Istream& is
) const
{
    while (!is.eof())
    {
        token t;
        readToken(t, is);
        if (t.isFunctionName())
        {
            if
            (
                t.functionNameToken() == ifEntry::typeName
             || t.functionNameToken() == ifeqEntry::typeName
            )
            {
                stack.append(filePos(is.name(), is.lineNumber()));
                skipUntil(stack, contextDict, endifEntry::typeName, is);
                stack.remove();
            }
            else if (t.functionNameToken() == endWord)
            {
                return;
            }
        }
    }

    FatalIOErrorInFunction(contextDict)
        << "Did not find matching " << endWord
        << " for " << typeName << " condition"
        << exit(FatalIOError);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::ifeqEntry::ifeqEntry
(
    const functionName& functionType,
    const dictionary& parentDict,
    const Istream& is,
    const tokenList& tokens
)
:
    functionEntry(functionType, parentDict, is, tokens)
{}


Foam::functionEntries::ifeqEntry::ifeqEntry
(
    const dictionary& parentDict,
    Istream& is
)
:
    ifeqEntry(typeName, parentDict, is, readArgList(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::ifeqEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    DynamicList<filePos> stack(10);
    return execute(stack, contextDict, contextDict, is);
}


bool Foam::functionEntries::ifeqEntry::execute
(
    const dictionary& contextDict,
    primitiveEntry& contextEntry,
    Istream& is
)
{
    DynamicList<filePos> stack(10);
    const ifeqEntry ifeqe(contextDict, is);

    return ifeqe.execute(stack, contextDict, contextEntry, is);
}


// ************************************************************************* //
