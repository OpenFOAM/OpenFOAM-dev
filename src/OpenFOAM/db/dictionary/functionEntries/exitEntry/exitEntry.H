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
    Foam::functionEntries::exitEntry

Description
    Exit the run without error

    providing the file name and line number of the #exit statement.
    \c \#exit does not require an argument.

See also
    Foam::functionEntries::printEntry

SourceFiles
    exitEntry.C

\*---------------------------------------------------------------------------*/

#ifndef exitEntry_H
#define exitEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                           Class exitEntry Declaration
\*---------------------------------------------------------------------------*/

class exitEntry
:
    public functionEntry
{

public:

    //- Runtime type information
    FunctionTypeName("#exit");


    // Constructors

        //- Construct from line number, dictionary and Istream
        exitEntry
        (
            const label lineNumber,
            const dictionary& parentDict,
            Istream& is
        );

        //- Copy construct
        exitEntry(const exitEntry&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new exitEntry(*this));
        }


    // Member Functions

        //- Expand the functionEntry into the contextDict
        virtual bool execute(dictionary& contextDict, Istream& is);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const exitEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
