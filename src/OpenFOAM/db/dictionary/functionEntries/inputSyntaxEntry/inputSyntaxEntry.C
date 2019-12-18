/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "inputSyntaxEntry.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionEntries::inputSyntaxEntry::typeName
(
    Foam::functionEntries::inputSyntaxEntry::typeName_()
);

// Don't lookup the debug switch here as the debug switch dictionary
// might include inputSyntax entries
int Foam::functionEntries::inputSyntaxEntry::debug(0);

// Read the default dictionary syntax from etc/controlDict if specified
Foam::functionEntries::inputSyntaxEntry::inputSyntax
    Foam::functionEntries::inputSyntaxEntry::defaultSyntax_
(
    Foam::debug::optimisationSwitches().found("inputSyntax")
  ? Foam::functionEntries::inputSyntaxEntry::syntax
    (
        Foam::debug::optimisationSwitches().lookup
        (
            "inputSyntax"
        )
    )
  : DOT
);

// Initialise the current dictionary syntax to the default
Foam::functionEntries::inputSyntaxEntry::inputSyntax
    Foam::functionEntries::inputSyntaxEntry::syntax_
(
    Foam::functionEntries::inputSyntaxEntry::defaultSyntax_
);


namespace Foam
{
namespace functionEntries
{
    addToMemberFunctionSelectionTable
    (
        functionEntry,
        inputSyntaxEntry,
        execute,
        dictionaryIstream
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::functionEntries::inputSyntaxEntry::inputSyntax
Foam::functionEntries::inputSyntaxEntry::syntax
(
    Istream& is
)
{
    word syntax(is);
    if (syntax == "slash")
    {
        return SLASH;
    }
    else if (syntax == "dot")
    {
        return DOT;
    }
    else
    {
        WarningInFunction
            << "unsupported input syntax'" << syntax
            << ", setting to default"
            << endl;

        return defaultSyntax_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::inputSyntaxEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    syntax_ = syntax(is);
    return true;
}


void Foam::functionEntries::inputSyntaxEntry::clear()
{
    syntax_ = defaultSyntax_;
}


bool Foam::functionEntries::inputSyntaxEntry::slash()
{
    return syntax_ == SLASH;
}


bool Foam::functionEntries::inputSyntaxEntry::dot()
{
    return syntax_ == DOT;
}


char Foam::functionEntries::inputSyntaxEntry::scopeChar()
{
    return syntax_ == SLASH ? '/' : '.';
}


// ************************************************************************* //
