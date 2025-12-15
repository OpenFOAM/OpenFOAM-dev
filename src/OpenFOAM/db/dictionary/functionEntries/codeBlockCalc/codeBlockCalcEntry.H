/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::functionEntries::codeBlockCalcEntry

Description

See also
    Foam::functionEntries::codeIncludeEntry
    Foam::functionEntries::codeStream

SourceFiles
    codeBlockCalcEntry.C

\*---------------------------------------------------------------------------*/

#ifndef codeBlockCalcEntry_H
#define codeBlockCalcEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

class codeBlockEntry;

/*---------------------------------------------------------------------------*\
                          Class codeBlockCalcEntry Declaration
\*---------------------------------------------------------------------------*/

class codeBlockCalcEntry
:
    public functionEntry
{
        //- Pointer to the codeBlockEntry which created this codeBlockCalcEntry
        const codeBlockEntry* codeBlockPtr_;

        //- Index of the function in the codeBlock
        //  corresponding to this codeBlockCalcEntry
        label index_;


    // Private Member Functions

        //- Perform the calculation and return the resulting string
        string calc(const dictionary& dict, Istream& is);


public:

    //- Runtime type information
    FunctionTypeName("#codeBlockCalc");


    // Constructors

        //- Construct from line number, dictionary and Istream
        codeBlockCalcEntry
        (
            const label lineNumber,
            const dictionary& parentDict,
            Istream& is
        );

        //- Copy construct
        codeBlockCalcEntry(const codeBlockCalcEntry&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new codeBlockCalcEntry(*this));
        }


    // Member Functions

        //- Expand the functionEntry into the contextDict
        virtual bool execute(dictionary& contextDict, Istream&);

        //- Expand the functionEntry into the contextEntry
        static bool execute
        (
            const dictionary& contextDict,
            primitiveEntry& contextEntry,
            Istream&
        );

        //- Dummy write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const codeBlockCalcEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
