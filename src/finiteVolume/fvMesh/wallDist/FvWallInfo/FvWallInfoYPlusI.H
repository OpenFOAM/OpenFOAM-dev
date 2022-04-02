/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "FvWallInfoYPlus.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class WallInfo, class Derived>
template<class TrackingData>
inline bool Foam::FvWallInfoYPlusBase<WallInfo, Derived>::update
(
    const point& pt,
    const FvWallInfoYPlusBase<WallInfo, Derived>& w2,
    const scalar tol,
    TrackingData& td
)
{
    FvWallInfoYPlusBase<WallInfo, Derived> copy(*this);

    bool result =
        FvWallInfoDataBase<WallInfo, scalar, Derived>::update(pt, w2, tol, td);

    if (result)
    {
        const scalar yPlus = sqrt(this->distSqr())/w2.data();

        if (yPlus >= td.yPlusCutOff)
        {
            *this = copy;
            result = false;
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class WallInfo, class Derived>
inline Foam::FvWallInfoYPlusBase<WallInfo, Derived>::FvWallInfoYPlusBase()
:
    FvWallInfoDataBase<WallInfo, scalar, Derived>()
{
    // Important: The value of yStar where the meshWave does not come.
    this->data() = 1.0;
}


// ************************************************************************* //
