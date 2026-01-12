/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::functionEntries::codeBlockDictEntry

Description

See also
    Foam::functionEntries::codeBlockStream

SourceFiles
    codeBlockDictEntry.C

\*---------------------------------------------------------------------------*/

#ifndef codeBlockDictEntry_H
#define codeBlockDictEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

class codeBlockEntry;

/*---------------------------------------------------------------------------*\
                          Class codeBlockDictEntry Declaration
\*---------------------------------------------------------------------------*/

class codeBlockDictEntry
:
    public functionEntry
{
        //- Pointer to the codeBlockEntry which created this codeBlockDictEntry
        const codeBlockEntry* codeBlockPtr_;

        //- Index of the function in the codeBlock
        //  corresponding to this codeBlockDictEntry
        label index_;


    // Private Member Functions

        //- Compile and execute the code to manipulate the contextDict
        bool resultStream(dictionary& contextDict, Istream& is);


public:

    //- Runtime type information
    FunctionTypeName("#codeBlockDict");


    // Constructors

        //- Construct from line number, dictionary and Istream
        codeBlockDictEntry
        (
            const label lineNumber,
            const dictionary& parentDict,
            Istream& is
        );

        //- Copy construct
        codeBlockDictEntry(const codeBlockDictEntry&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new codeBlockDictEntry(*this));
        }


    // Member Functions

        //- Expand the functionEntry into the contextDict
        virtual bool execute(dictionary& contextDict, Istream&);

        //- Dummy write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const codeBlockDictEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
