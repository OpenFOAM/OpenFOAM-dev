/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2026 OpenFOAM Foundation
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

#include "unitSet.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline Foam::unitSet Foam::unitSet::newUnitless()
{
    return unitSet(dimensionSet(0, 0, 0, 0, 0), 0, 0, 1);
}


inline Foam::unitSet Foam::unitSet::newAny()
{
    return unitSet(dimensionSet(0, 0, 0, 0, 0), 0, 0, 0);
}


inline Foam::unitSet Foam::unitSet::newNone()
{
    return unitSet(dimensionSet(0, 0, 0, 0, 0), 0, 0, -1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::dimensionSet& Foam::unitSet::dimensions() const
{
    return dimensions_;
}


inline bool Foam::unitSet::dimensionless() const
{
    return dimensions_.dimensionless();
}


inline bool Foam::unitSet::unitless() const
{
    return compare(*this, newUnitless(), false);
}


inline bool Foam::unitSet::any() const
{
    return multiplier_ == 0;
}


inline bool Foam::unitSet::none() const
{
    return multiplier_ == -1;
}


inline bool Foam::unitSet::standard() const
{
    return multiplier_ == 1 || any() || none();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::scalar Foam::unitSet::operator[]
(
    const dimensionSet::dimensionType type
) const
{
    return dimensions_[type];
}


inline Foam::scalar Foam::unitSet::operator[]
(
    const dimlessUnitType type
) const
{
    return exponents_[static_cast<int>(type)];
}


// ************************************************************************* //
