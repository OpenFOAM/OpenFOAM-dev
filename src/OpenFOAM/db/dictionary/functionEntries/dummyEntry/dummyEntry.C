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

\*---------------------------------------------------------------------------*/

#include "dummyEntry.H"
#include "includeFuncEntry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeName(dummyEntry);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::dummyEntry::dummyEntry
(
    const functionName& functionType,
    const dictionary& parentDict,
    Istream& is
)
:
    functionEntry(functionType, is.lineNumber(), parentDict)
{
    if (functionType == includeFuncEntry::typeName_())
    {
        cerr<< "--> FOAM Error: Found " << functionType
            << " while reading system/controlDict"
            << std::endl
            << "    Move the functions entries into the system/functions file"
            << std::endl;
        std::exit(1);
    }

    // Get the rest of the line and discard
    string line;
    dynamic_cast<ISstream&>(is).getLine(line);
}


// ************************************************************************* //
