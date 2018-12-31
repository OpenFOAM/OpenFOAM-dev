/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
#include "stringOps.H"
#include "ifEntry.H"
#include "Switch.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(ifeqEntry, 0);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        ifeqEntry,
        execute,
        dictionaryIstream,
        ifeq
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionEntries::ifeqEntry::readToken(token& t, Istream& is)
{
    // Skip dummy tokens - avoids entry::getKeyword consuming #else, #endif
    do
    {
        if
        (
            is.read(t).bad()
         || is.eof()
         || !t.good()
        )
        {
            return;
        }
    }
    while (t == token::END_STATEMENT);
}


Foam::token Foam::functionEntries::ifeqEntry::expand
(
    const dictionary& dict,
    const string& keyword,
    const token& t
)
{
    if (keyword[0] == '$')
    {
        word varName = keyword(1, keyword.size()-1);

        // lookup the variable name in the given dictionary
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
            string expanded(keyword);
            stringOps::inplaceExpand(expanded, dict, true, true);

            // Re-form as a string token so we can compare to string
            return token(expanded, t.lineNumber());
        }
    }
    else if (!t.isString())
    {
        // Re-form as a string token so we can compare to string
        return token(keyword, t.lineNumber());
    }
    else
    {
        return t;
    }
}


Foam::token Foam::functionEntries::ifeqEntry::expand
(
    const dictionary& dict,
    const token& t
)
{
    if (t.isWord())
    {
        return expand(dict, t.wordToken(), t);
    }
    else if (t.isVariable())
    {
        return expand(dict, t.stringToken(), t);
    }
    else if (t.isString())
    {
        return expand(dict, t.stringToken(), t);
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
)
{
    const bool eqType = (t1.type() == t2.type());

    switch (t1.type())
    {
        case token::UNDEFINED:
            return eqType;

        case token::PUNCTUATION:
            return (eqType && t1.pToken() == t2.pToken());

        case token::WORD:
            if (eqType)
            {
                return t1.wordToken() == t2.wordToken();
            }
            else if (t2.isString())
            {
                return t1.wordToken() == t2.stringToken();
            }
            else
            {
                return false;
            }

        case token::STRING:
        case token::VARIABLE:
        case token::VERBATIMSTRING:
            if (eqType)
            {
                return t1.stringToken() == t2.stringToken();
            }
            else if (t2.isWord())
            {
                return t1.stringToken() == t2.wordToken();
            }
            else
            {
                return false;
            }

        case token::LABEL:
            if (eqType)
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
    const dictionary& parentDict,
    const word& endWord,
    Istream& is
)
{
    while (!is.eof())
    {
        token t;
        readToken(t, is);
        if (t.isWord())
        {
            if (t.wordToken() == "#if" || t.wordToken() == "#ifeq")
            {
                stack.append(filePos(is.name(), is.lineNumber()));
                skipUntil(stack, parentDict, "#endif", is);
                stack.remove();
            }
            else if (t.wordToken() == endWord)
            {
                return;
            }
        }
    }

    FatalIOErrorInFunction(parentDict)
        << "Did not find matching " << endWord << exit(FatalIOError);
}


bool Foam::functionEntries::ifeqEntry::evaluate
(
    const bool doIf,
    DynamicList<filePos>& stack,
    dictionary& parentDict,
    Istream& is
)
{
    while (!is.eof())
    {
        token t;
        readToken(t, is);

        if (t.isWord() && t.wordToken() == "#ifeq")
        {
            // Recurse to evaluate
            execute(stack, parentDict, is);
        }
        else if (t.isWord() && t.wordToken() == "#if")
        {
            // Recurse to evaluate
            ifEntry::execute(stack, parentDict, is);
        }
        else if
        (
            doIf
         && t.isWord()
         && (t.wordToken() == "#else" || t.wordToken() == "#elif")
        )
        {
            // Now skip until #endif
            skipUntil(stack, parentDict, "#endif", is);
            stack.remove();
            break;
        }
        else if (t.isWord() && t.wordToken() == "#endif")
        {
            stack.remove();
            break;
        }
        else
        {
            is.putBack(t);
            bool ok = entry::New(parentDict, is);
            if (!ok)
            {
                return false;
            }
        }
    }
    return true;
}


bool Foam::functionEntries::ifeqEntry::execute
(
    const bool doIf,
    DynamicList<filePos>& stack,
    dictionary& parentDict,
    Istream& is
)
{
    if (doIf)
    {
        evaluate(true, stack, parentDict, is);
    }
    else
    {
        // Fast-forward to #else
        token t;
        while (!is.eof())
        {
            readToken(t, is);
            if
            (
                t.isWord()
             && (t.wordToken() == "#if" || t.wordToken() == "#ifeq")
            )
            {
                stack.append(filePos(is.name(), is.lineNumber()));
                skipUntil(stack, parentDict, "#endif", is);
                stack.remove();
            }
            else if (t.isWord() && t.wordToken() == "#else")
            {
                break;
            }
            else if (t.isWord() && t.wordToken() == "#elif")
            {
                // const label lineNo = is.lineNumber();

                // Read line
                string line;
                dynamic_cast<ISstream&>(is).getLine(line);
                line += ';';
                IStringStream lineStream(line);
                const primitiveEntry e("ifEntry", parentDict, lineStream);
                const Switch doIf(e.stream());

                if (doIf)
                {
                    // Info<< "Using #elif " << doIf << " at line " << lineNo
                    //     << " in file " << is.name() << endl;
                    break;
                }
            }
            else if (t.isWord() && t.wordToken() == "#endif")
            {
                stack.remove();
                break;
            }
        }

        if (t.wordToken() == "#else")
        {
            // Evaluate until we hit #endif
            evaluate(false, stack, parentDict, is);
        }
        else if (t.wordToken() == "#elif")
        {
            // Evaluate until we hit #else or #endif
            evaluate(true, stack, parentDict, is);
        }
    }
    return true;
}


bool Foam::functionEntries::ifeqEntry::execute
(
    DynamicList<filePos>& stack,
    dictionary& parentDict,
    Istream& is
)
{
    const label nNested = stack.size();

    stack.append(filePos(is.name(), is.lineNumber()));

    // Read first token and expand any string
    token cond1(is);
    cond1 = expand(parentDict, cond1);

    // Read second token and expand any string
    token cond2(is);
    cond2 = expand(parentDict, cond2);

    const bool equal = equalToken(cond1, cond2);

    // Info<< "Using #" << typeName << " " << cond1
    //     << " == " << cond2
    //     << " at line " << stack.last().second()
    //     << " in file " <<  stack.last().first() << endl;

    bool ok = ifeqEntry::execute(equal, stack, parentDict, is);

    if (stack.size() != nNested)
    {
        FatalIOErrorInFunction(parentDict)
            << "Did not find matching #endif for condition starting"
            << " at line " << stack.last().second()
            << " in file " <<  stack.last().first() << exit(FatalIOError);
    }

    return ok;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::ifeqEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    DynamicList<filePos> stack(10);
    return execute(stack, parentDict, is);
}


// ************************************************************************* //
