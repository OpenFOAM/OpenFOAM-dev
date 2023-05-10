/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::scalar Foam::integerPow(const scalar x, const label e)
{
    if (e == 0)
    {
        return 1;
    }

    scalar xx = e > 0 ? x : 1/x;

    scalar y = 1;

    for (label i = e; i != 0; i /= 2)
    {
        if (i % 2)
        {
            y *= xx;
        }

        xx *= xx;
    }

    return y;
}


inline Foam::scalar Foam::integerRoot(const scalar x, const label e)
{
    switch (e)
    {
        case -4: return 1/pow025(x);
        case -3: return 1/cbrt(x);
        case -2: return 1/sqrt(x);
        case -1: return 1/x;
        case 0: return NaN;
        case 1: return x;
        case 2: return sqrt(x);
        case 3: return cbrt(x);
        case 4: return pow025(x);
    }

    return pow(x, scalar(1)/e);
}


// ************************************************************************* //
