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

#include "ifEntry.H"
#include "elifEntry.H"
#include "elseEntry.H"
#include "endifEntry.H"
#include "Switch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Context>
bool Foam::functionEntries::ifeqEntry::evaluate
(
    const bool doIf,
    DynamicList<filePos>& stack,
    const dictionary& contextDict,
    Context& context,
    Istream& is
) const
{
    while (!is.eof())
    {
        token t;
        readToken(t, is);

        if (t.isFunctionName() && t.functionNameToken() == ifeqEntry::typeName)
        {
            // Recurse to evaluate
            execute(stack, contextDict, context, is);
        }
        else if
        (
            t.isFunctionName()
         && t.functionNameToken() == ifEntry::typeName
        )
        {
            // Recurse to evaluate
            const ifEntry ife(contextDict, is);
            ife.execute(stack, contextDict, context, is);
        }
        else if
        (
            doIf
         && t.isFunctionName()
         && (
                t.functionNameToken() == elseEntry::typeName
             || t.functionNameToken() == elifEntry::typeName
            )
        )
        {
            // Now skip until #endif
            skipUntil(stack, contextDict, endifEntry::typeName, is);
            stack.remove();
            break;
        }
        else if
        (
            t.isFunctionName()
         && t.functionNameToken() == endifEntry::typeName
        )
        {
            stack.remove();
            break;
        }
        else
        {
            const bool ok = insert(contextDict, context, t, is);
            if (!ok)
            {
                return false;
            }
        }
    }

    return true;
}


template<class Context>
bool Foam::functionEntries::ifeqEntry::execute
(
    const bool doIf,
    DynamicList<filePos>& stack,
    const dictionary& contextDict,
    Context& context,
    Istream& is
) const
{
    if (doIf)
    {
        evaluate(true, stack, contextDict, context, is);
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
                    t.functionNameToken() == ifEntry::typeName
                 || t.functionNameToken() == ifeqEntry::typeName
                )
                {
                    stack.append(filePos(is.name(), is.lineNumber()));
                    skipUntil(stack, contextDict, endifEntry::typeName, is);
                    stack.remove();
                }
                else if (t.functionNameToken() == elseEntry::typeName)
                {
                    break;
                }
                else if (t.functionNameToken() == elifEntry::typeName)
                {
                    // Read elif argument by constructing an ifEntry
                    const ifEntry ife(contextDict, is);
                    const string arg
                    (
                        ife[1].stringToken() + char(token::END_STATEMENT)
                    );
                    IStringStream argStream(arg);
                    argStream.lineNumber() = ife[1].lineNumber();
                    const primitiveEntry e("elifEntry", contextDict, argStream);
                    const Switch doIf(e.stream());

                    if (doIf)
                    {
                        break;
                    }
                }
                else if (t.functionNameToken() == endifEntry::typeName)
                {
                    stack.remove();
                    break;
                }
            }
        }

        if (t.functionNameToken() == elseEntry::typeName)
        {
            // Evaluate until we hit #endif
            evaluate(false, stack, contextDict, context, is);
        }
        else if (t.functionNameToken() == elifEntry::typeName)
        {
            // Evaluate until we hit #else or #endif
            evaluate(true, stack, contextDict, context, is);
        }
    }
    return true;
}


template<class Context>
bool Foam::functionEntries::ifeqEntry::execute
(
    DynamicList<filePos>& stack,
    const dictionary& contextDict,
    Context& context,
    Istream& is
) const
{
    const label nNested = stack.size();

    stack.append(filePos(is.name(), is.lineNumber()));

    // Read first token and expand if a variable
    token cond1(operator[](1));
    cond1 = expand(contextDict, cond1);

    // Read second token and expand if a variable
    token cond2(operator[](2));
    cond2 = expand(contextDict, cond2);

    const bool equal = equalToken(cond1, cond2);

    // Info<< "Using #" << typeName << " " << cond1
    //     << " == " << cond2
    //     << " at line " << stack.last().second()
    //     << " in file " <<  stack.last().first() << endl;

    bool ok = ifeqEntry::execute(equal, stack, contextDict, context, is);

    if (stack.size() != nNested)
    {
        FatalIOErrorInFunction(contextDict)
            << "Did not find matching #endif for "
            << typeName << " condition starting"
            << " at line " << stack.last().second()
            << " in file " <<  stack.last().first() << exit(FatalIOError);
    }

    return ok;
}


// ************************************************************************* //
