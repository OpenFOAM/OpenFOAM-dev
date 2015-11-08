/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "error.H"
#include "uLabel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if WM_LABEL_SIZE == 32
const char* const Foam::pTraits<uint64_t>::typeName = "uint64";
const char* const Foam::pTraits<uint32_t>::typeName = "uLabel";
#elif WM_LABEL_SIZE == 64
const char* const Foam::pTraits<uint64_t>::typeName = "uLabel";
const char* const Foam::pTraits<uint32_t>::typeName = "uint32";
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::uLabel Foam::pow(uLabel a, uLabel b)
{
    uLabel ans = 1;
    for (uLabel i=0; i<b; i++)
    {
        ans *= a;
    }

    #ifdef FULLDEBUG
    if (b < 0)
    {
        FatalErrorInFunction
            << "negative value for b is not supported"
            << abort(FatalError);
    }
    #endif

    return ans;
}


Foam::uLabel Foam::factorial(uLabel n)
{
    static uLabel factTable[13] =
    {
        1, 1, 2, 6, 24, 120, 720, 5040, 40320,
        362880, 3628800, 39916800, 479001600
    };

    #ifdef FULLDEBUG
    if (n > 12 && n < 0)
    {
        FatalErrorInFunction
            << "n value out of range"
            << abort(FatalError);
    }
    #endif

    return factTable[n];
}


// ************************************************************************* //
