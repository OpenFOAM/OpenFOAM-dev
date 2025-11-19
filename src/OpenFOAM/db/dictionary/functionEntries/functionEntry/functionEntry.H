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

Namespace
    Foam::functionEntries

Description
    Namespace for containing a functionEntry.


Class
    Foam::functionEntry

Description
    A functionEntry causes entries to be added/manipulated on the specified
    dictionary given an input stream.

    In dictionaries, a \c '\#' sigil is typically used for a functionEntry.

SourceFiles
    functionEntry.C
    functionEntryIO.C

\*---------------------------------------------------------------------------*/

#ifndef functionEntry_H
#define functionEntry_H

#include "primitiveEntry.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

/*---------------------------------------------------------------------------*\
                        Class functionEntry Declaration
\*---------------------------------------------------------------------------*/

class functionEntry
:
    public primitiveEntry
{
    // Private Member Functions

        //- Read the rest of the argument list following the part attached
        // to the function name
        void readRestOfArgs(string& fNameArgs, Istream& is) const;


protected:

    // Constructors

        //- Construct from function type and parent dictionary
        functionEntry
        (
            const functionName&,
            const label lineNumber,
            const dictionary&
        );

        //- Construct from function type, parent dictionary and a token
        functionEntry
        (
            const functionName&,
            const label lineNumber,
            const dictionary&,
            const Istream& is,
            const token&
        );

        //- Construct from function type, parent dictionary and tokenList
        functionEntry
        (
            const functionName&,
            const label lineNumber,
            const dictionary&,
            const Istream& is,
            const tokenList&
        );

        //- Copy construct
        functionEntry(const functionEntry&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const = 0;


    // Protected member functions

        tokenList readArgList
        (
            const functionName& functionType,
            Istream& is,
            const bool optional = false
        ) const;

        //- Read the function name and optional argument list
        //  The argument list can be on one or more lines.
        tokenList readFuncNameArgList
        (
            const functionName& functionType,
            Istream& is
        ) const;

        //- Read the filen name and optional argument list
        //  The argument list can be on one or more lines.
        tokenList readFileNameArgList
        (
            const functionName& functionType,
            Istream& is
        ) const;

        //- Read entry starting from token t and insert into contextDict
        static bool insert
        (
            const dictionary& parentDict,
            dictionary& contextDict,
            const token& t,
            Istream& is
        );

        //- Insert token t into contextEntry, expanding as required
        //  which may read subsequent tokens from Istream
        static bool insert
        (
            const dictionary& parentDict,
            primitiveEntry& contextEntry,
            const token& t,
            Istream& is
        );

        //- Read entry from the given string and insert into contextDict
        static bool insert
        (
            dictionary& contextDict,
            const string&
        );

        //- Read tokens from the given string and insert into contextEntry
        static bool insert
        (
            const dictionary& parentDict,
            primitiveEntry& contextEntry,
            const string&
        );


public:

    FunctionTypeName("functionEntry");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            functionEntry,
            dictionary,
            (
                const label lineNumber,
                const dictionary& parentDict,
                Istream& is
            ),
            (lineNumber, parentDict, is)
        );


    // Selectors

        //- Construct from Istream and insert into dictionary
        static autoPtr<functionEntry> New
        (
            const keyType& functionName,
            const dictionary& parentDict,
            Istream& is
        );


    // Member Function Selectors

        //- Expand the functionEntry into the contextDict
        virtual bool execute(dictionary& contextDict, Istream& is) = 0;

        declareMemberFunctionSelectionTable
        (
            bool,
            functionEntry,
            execute,
            primitiveEntryIstream,
            (
                const dictionary& contextDict,
                primitiveEntry& contextEntry,
                Istream& is
            ),
            (contextDict, entry, is)
        );

        //- Expand the functionEntry into the contextEntry
        static bool execute
        (
            const word& functionName,
            const dictionary& contextDict,
            primitiveEntry& contextEntry,
            Istream&
        );

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const functionEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
