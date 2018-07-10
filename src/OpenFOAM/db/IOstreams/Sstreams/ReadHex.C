/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Read a non-delimited hex label

\*---------------------------------------------------------------------------*/

#include "ReadHex.H"
#include <cctype>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
T Foam::ReadHex(ISstream& is)
{
    // Takes into account that 'a' (or 'A') is 10
    static const int alphaOffset = toupper('A') - 10;
    // Takes into account that '0' is 0
    static const int zeroOffset = int('0');

    char c = 0;

    // Get next non-whitespace character
    while (is.get(c) && isspace(c))
    {}

    T result = 0;
    do
    {
        if (isspace(c) || c == 0) break;

        if (!isxdigit(c))
        {
            FatalIOErrorInFunction(is)
                << "Illegal hex digit: '" << c << "'"
                << exit(FatalIOError);
        }

        result <<= 4;

        if (isdigit(c))
        {
            result += int(c) - zeroOffset;
        }
        else
        {
            result += toupper(c) - alphaOffset;
        }
    } while (is.get(c));

    return result;
}


// ************************************************************************* //
