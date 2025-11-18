/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "includeFvConstraintEntry.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(includeFvConstraintEntry, 0);

    addToRunTimeSelectionTable
    (
        functionEntry,
        includeFvConstraintEntry,
        dictionary
    );
}
}


Foam::fileName
Foam::functionEntries::includeFvConstraintEntry::fvConstraintDictPath
(
    "caseDicts/fvConstraints"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::includeFvConstraintEntry::includeFvConstraintEntry
(
    const dictionary& parentDict,
    Istream& is
)
:
    includeFuncEntry(typeName, parentDict, is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeFvConstraintEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    return readConfigFile
    (
        "constraint",
        funcNameArgs(),
        contextDict,
        fvConstraintDictPath,
        "system"
    );
}


// ************************************************************************* //
