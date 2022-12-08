/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "word.H"
#include "memberFunctionSelectionTables.H"
#include "primitiveEntry.H"

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

        //- Read line as word token
        static token readLine(Istream& is);


protected:

        //- Read the function name and optional argument list
        //  and return as a string.
        //  The optional argument list can be on one or more lines.
        static string readFuncNameArgs(Istream& is);

        //- Insert the given string in a sub-dict context
        static bool insert
        (
            dictionary& parentDict,
            const string&
        );

        //- Insert the given string in a primitiveEntry context
        static bool insert
        (
            const dictionary& parentDict,
            primitiveEntry& thisEntry,
            const string&
        );


public:

    // Constructors

        //- Construct from keyword, parent dictionary and Istream
        functionEntry(const word&, const dictionary&, Istream&);

        //- Copy construct
        functionEntry(const functionEntry&) = default;

        //- Clone
        autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new functionEntry(*this));
        }


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            bool,
            functionEntry,
            execute,
            dictionaryIstream,
            (
                dictionary& parentDict,
                Istream& is
            ),
            (parentDict, is)
        );

        //- Execute the functionEntry in a sub-dict context
        static bool execute
        (
            const word& functionName,
            dictionary& parentDict,
            Istream&
        );

        declareMemberFunctionSelectionTable
        (
            bool,
            functionEntry,
            execute,
            primitiveEntryIstream,
            (
                const dictionary& parentDict,
                primitiveEntry& entry,
                Istream& is
            ),
            (parentDict, entry, is)
        );

        //- Execute the functionEntry in a primitiveEntry context
        static bool execute
        (
            const word& functionName,
            const dictionary& parentDict,
            primitiveEntry&,
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
