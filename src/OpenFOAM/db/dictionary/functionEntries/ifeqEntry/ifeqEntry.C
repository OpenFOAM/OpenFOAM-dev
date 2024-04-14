/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        ifeqEntry,
        execute,
        dictionaryIstream
    );
}
}

const Foam::functionName Foam::functionEntries::ifeqEntry::ifName("#if");
const Foam::functionName Foam::functionEntries::ifeqEntry::ifeqName("#ifeq");
const Foam::functionName Foam::functionEntries::ifeqEntry::elifName("#elif");
const Foam::functionName Foam::functionEntries::ifeqEntry::elseName("#else");
const Foam::functionName Foam::functionEntries::ifeqEntry::endifName("#endif");

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
    const token& t
)
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
    const dictionary& parentDict,
    const functionName& endWord,
    Istream& is
)
{
    while (!is.eof())
    {
        token t;
        readToken(t, is);
        if (t.isFunctionName())
        {
            if
            (
                t.functionNameToken() == ifName
             || t.functionNameToken() == ifeqName
            )
            {
                stack.append(filePos(is.name(), is.lineNumber()));
                skipUntil(stack, parentDict, endifName, is);
                stack.remove();
            }
            else if (t.functionNameToken() == endWord)
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

        if (t.isFunctionName() && t.functionNameToken() == ifeqName)
        {
            // Recurse to evaluate
            execute(stack, parentDict, is);
        }
        else if (t.isFunctionName() && t.functionNameToken() == ifName)
        {
            // Recurse to evaluate
            ifEntry::execute(stack, parentDict, is);
        }
        else if
        (
            doIf
         && t.isFunctionName()
         && (
                t.functionNameToken() == elseName
             || t.functionNameToken() == elifName
            )
        )
        {
            // Now skip until #endif
            skipUntil(stack, parentDict, endifName, is);
            stack.remove();
            break;
        }
        else if (t.isFunctionName() && t.functionNameToken() == endifName)
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

            if (t.isFunctionName())
            {
                if
                (
                    t.functionNameToken() == ifName
                 || t.functionNameToken() == ifeqName
                )
                {
                    stack.append(filePos(is.name(), is.lineNumber()));
                    skipUntil(stack, parentDict, endifName, is);
                    stack.remove();
                }
                else if (t.functionNameToken() == elseName)
                {
                    break;
                }
                else if (t.functionNameToken() == elifName)
                {
                    // Read line
                    string line;
                    dynamic_cast<ISstream&>(is).getLine(line);
                    line += ';';
                    IStringStream lineStream(line);
                    const primitiveEntry e("ifEntry", parentDict, lineStream);
                    const Switch doIf(e.stream());

                    if (doIf)
                    {
                        // Info<< "Using #elif " << doIf
                        //     << " at line " << lineNo
                        //     << " in file " << is.name() << endl;
                        break;
                    }
                }
                else if (t.functionNameToken() == endifName)
                {
                    stack.remove();
                    break;
                }
            }
        }

        if (t.functionNameToken() == elseName)
        {
            // Evaluate until we hit #endif
            evaluate(false, stack, parentDict, is);
        }
        else if (t.functionNameToken() == elifName)
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

    // Read first token and expand if a variable
    token cond1(is);
    cond1 = expand(parentDict, cond1);

    // Read second token and expand if a variable
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
