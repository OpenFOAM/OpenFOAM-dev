/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "HashTable.H"
#include "uLabel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(HashTableCore, 0);
}

const Foam::label Foam::HashTableCore::maxTableSize
(
    Foam::HashTableCore::canonicalSize
    (
        Foam::labelMax/2
    )
);


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::HashTableCore::canonicalSize(const label size)
{
    if (size < 1)
    {
        return 0;
    }

    // enforce power of two
    uLabel goodSize = size;

    if (goodSize & (goodSize - 1))
    {
        // brute-force is fast enough
        goodSize = 1;
        while (goodSize < unsigned(size))
        {
            goodSize <<= 1;
        }
    }

    return goodSize;
}


// ************************************************************************* //
