/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class primitiveEntry;
class dictionary;

/*---------------------------------------------------------------------------*\
                           Class functionEntry Declaration
\*---------------------------------------------------------------------------*/

class functionEntry
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionEntry(const functionEntry&);

        //- Disallow default bitwise assignment
        void operator=(const functionEntry&);


public:

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


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
