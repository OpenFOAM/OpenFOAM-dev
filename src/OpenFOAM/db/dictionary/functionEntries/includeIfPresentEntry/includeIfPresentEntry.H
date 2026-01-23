/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::functionEntries::includeIfPresentEntry

Description
    Specify a file to include if it exists. Expects a single string to follow.

    The \c \#includeIfPresent directive is similar to the \c \#include
    directive, but does not generate an error if the file does not exist.

    Additionally optional arguments may be specified which are available for
    substitution in the included file e.g.
    \verbatim
        #includeIfPresent "regionFunctions"(region = heatedWall)
    \endverbatim

See also
    Foam::functionEntries::includeEntry

SourceFiles
    includeIfPresentEntry.C

\*---------------------------------------------------------------------------*/

#ifndef includeIfPresentEntry_H
#define includeIfPresentEntry_H

#include "includeEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                    Class includeIfPresentEntry Declaration
\*---------------------------------------------------------------------------*/

class includeIfPresentEntry
:
    public includeEntry
{

public:

    //- Runtime type information
    FunctionTypeName("#includeIfPresent");


    // Constructors

        //- Construct from line number, dictionary and Istream
        includeIfPresentEntry
        (
            const label lineNumber,
            const dictionary& parentDict,
            Istream& is
        );

        //- Copy construct
        includeIfPresentEntry(const includeIfPresentEntry&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new includeIfPresentEntry(*this));
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


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const includeIfPresentEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
