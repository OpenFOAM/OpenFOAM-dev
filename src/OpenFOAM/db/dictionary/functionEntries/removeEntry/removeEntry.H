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
#include "wordReList.H"

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
    // Private Member Data

        //- The list of patterns of entries to remove
        wordReList patterns_;


    // Private Member Functions

        //- Read a single argument or list
        tokenList readArgOrList
        (
            const functionName& functionType,
            Istream& is
        ) const;


public:

    //- Runtime type information
    FunctionTypeName("#remove");


    // Constructors

        //- Construct from dictionary and Istream
        removeEntry(const dictionary& parentDict, Istream& is);

        //- Copy construct
        removeEntry(const removeEntry&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new removeEntry(*this));
        }


    // Member Functions

        //- Expand the functionEntry into the contextDict
        virtual bool execute(dictionary& contextDict, Istream&);


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
