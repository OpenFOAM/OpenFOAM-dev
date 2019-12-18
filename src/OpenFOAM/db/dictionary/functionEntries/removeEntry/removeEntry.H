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

Class
    Foam::functionEntries::removeEntry

Description
    Remove a dictionary entry.

    The \c \#remove directive takes a list or a single wordRe.
    For example,
    \verbatim
        #remove entry0
        #remove ( entry1 entry2 entry3 otherEntry )
        #remove "entry[1-3]"
        #remove ( "entry[1-3]" otherEntry )
    \endverbatim

    The removal only occurs in the current context.
    Removing sub-entries or parent entries is not supported.

SourceFiles
    removeEntry.C

\*---------------------------------------------------------------------------*/

#ifndef removeEntry_H
#define removeEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                           Class removeEntry Declaration
\*---------------------------------------------------------------------------*/

class removeEntry
:
    public functionEntry
{

public:

    //- Runtime type information
    ClassName("remove");


    // Constructors

        //- Disallow default bitwise copy construction
        removeEntry(const removeEntry&) = delete;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const removeEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
