/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "barycentric2D.H"
#include "Random.H"
#include "cachedRandom.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::barycentric2D barycentric2D01
(
    Foam::scalar s,
    Foam::scalar t
)
{
    // Transform the random point in the unit square to a random point in the
    // unit tri by reflecting across the diagonal

    if (s + t > 1)
    {
        s = 1 - s;
        t = 1 - t;
    }

    return Foam::barycentric2D(1 - s - t, s, t);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::barycentric2D Foam::barycentric2D01(Random& rndGen)
{
    return
        ::barycentric2D01
        (
            rndGen.scalar01(),
            rndGen.scalar01()
        );
}


Foam::barycentric2D Foam::barycentric2D01(cachedRandom& rndGen)
{
    return
        ::barycentric2D01
        (
            rndGen.sample01<scalar>(),
            rndGen.sample01<scalar>()
        );
}


// ************************************************************************* //
