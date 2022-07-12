/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

Description
    Class for handling debugging switches.

\*---------------------------------------------------------------------------*/

#include "runTimeSelectionToC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

bool Foam::debug::enableRunTimeSelectionToC = false;

Foam::HashTable<Foam::Tuple2<Foam::word, Foam::wordHashSet>>
    Foam::debug::runTimeSelectionToC;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::addToRunTimeSelectionTableToC
(
    const word& baseType,
    const word& baseTypeName,
    const word& thisTypeName
)
{
    if (debug::enableRunTimeSelectionToC)
    {
        if (!debug::runTimeSelectionToC.found(baseType))
        {
            debug::runTimeSelectionToC.insert
            (
                baseType,
                Tuple2<word, wordHashSet>(baseTypeName, wordHashSet())
            );
        }

        debug::runTimeSelectionToC[baseType].second().insert(thisTypeName);
    }

    return debug::enableRunTimeSelectionToC;
}


// ************************************************************************* //
