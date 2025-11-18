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

#include "includeFvModelEntry.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(includeFvModelEntry, 0);
    addToRunTimeSelectionTable(functionEntry, includeFvModelEntry, dictionary);
}
}


Foam::fileName Foam::functionEntries::includeFvModelEntry::fvModelDictPath
(
    "caseDicts/fvModels"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::includeFvModelEntry::includeFvModelEntry
(
    const dictionary& parentDict,
    Istream& is
)
:
    includeFuncEntry(typeName, parentDict, is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeFvModelEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    return readConfigFile
    (
        "model",
        funcNameArgs(),
        contextDict,
        fvModelDictPath,
        "constant"
    );
}


// ************************************************************************* //
