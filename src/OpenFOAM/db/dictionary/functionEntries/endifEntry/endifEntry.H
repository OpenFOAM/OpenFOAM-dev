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
    Foam::functionEntries::endifEntry

Description
    Part of the #if...#endif clause

See also
    Foam::functionEntries::ifEntry
    Foam::functionEntries::ifeqEntry

SourceFiles
    endifEntry.C

\*---------------------------------------------------------------------------*/

#ifndef endifEntry_H
#define endifEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                           Class endifEntry Declaration
\*---------------------------------------------------------------------------*/

class endifEntry
:
    public functionEntry
{

public:

    //- Runtime type information
    FunctionTypeName("#endif");


    // Constructors

        //- Construct from dictionary and Istream
        endifEntry(const dictionary& parentDict, Istream& is);

        //- Disallow default bitwise copy construction
        endifEntry(const endifEntry&) = delete;


    // Member Functions


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const endifEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
