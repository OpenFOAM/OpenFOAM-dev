/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "int32.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const int32_t Foam::pTraits<int32_t>::zero = 0;
const int32_t Foam::pTraits<int32_t>::one = 1;
const int32_t Foam::pTraits<int32_t>::min = INT32_MIN;
const int32_t Foam::pTraits<int32_t>::max = INT32_MAX;
const int32_t Foam::pTraits<int32_t>::rootMin = pTraits<int32_t>::min;
const int32_t Foam::pTraits<int32_t>::rootMax = pTraits<int32_t>::max;

const char* const Foam::pTraits<int32_t>::componentNames[] = { "" };

Foam::pTraits<int32_t>::pTraits(const int32_t& p)
:
    p_(p)
{}

Foam::pTraits<int32_t>::pTraits(Istream& is)
{
    is >> p_;
}


// ************************************************************************* //
