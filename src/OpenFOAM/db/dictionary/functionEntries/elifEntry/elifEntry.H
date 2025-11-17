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
    Foam::functionEntries::elifEntry

Description
    Conditional parsing of and within dictionary entries.

    Part of the #if...#elif...#endif

    E.g.
    \verbatim
        Uinlet 15;

        #if (#calc "${Uinlet} < 10")
            ..
        #elif (#calc "${Uinlet} > 20")
            ..
        #endif

        pOutlet
            #if (#calc "${Uinlet} < 10") 0
            #elif (#calc "${Uinlet} > 20") 100
            #else 10
            #endif;
    \endverbatim

    The condition should return any boolean representation supported by
    \c Switch i.e 0, 1, true, false, on, off, yes, no.

See also
    Foam::functionEntries::ifEntry
    Foam::functionEntries::ifeqEntry
    Foam::functionEntries::elseEntry
    Foam::functionEntries::endifEntry
    Foam::Switch

SourceFiles
    elifEntry.C

\*---------------------------------------------------------------------------*/

#ifndef elifEntry_H
#define elifEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                           Class elifEntry Declaration
\*---------------------------------------------------------------------------*/

class elifEntry
:
    public functionEntry
{

public:

    //- Runtime type information
    FunctionTypeName("#elif");


    // Constructors

        //- Construct from dictionary and Istream
        elifEntry(const dictionary& parentDict, Istream& is);

        //- Disallow default bitwise copy construction
        elifEntry(const elifEntry&) = delete;


    // Member Functions


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const elifEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
