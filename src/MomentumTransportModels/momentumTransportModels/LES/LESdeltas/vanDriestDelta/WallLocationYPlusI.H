/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "WallLocationYPlus.H"

// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

template<class WallLocation>
template<class TrackingYPlus>
inline bool Foam::WallLocationYPlus<WallLocation>::update
(
    const point& pt,
    const WallLocationYPlus<WallLocation>& w2,
    const scalar tol,
    TrackingYPlus& td
)
{
    WallLocationYPlus<WallLocation> copy(*this);

    bool result =
        WallLocationData<WallLocation, scalar>::update(pt, w2, tol, td);

    if (result && sqrt(this->distSqr())/w2.data() >= td.yPlusCutOff)
    {
        *this = copy;
        result = false;
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class WallLocation>
inline Foam::WallLocationYPlus<WallLocation>::WallLocationYPlus()
:
    WallLocationData<WallLocation, scalar>()
{
    // Important: The value of yStar where the meshWave does not come.
    this->data() = 1.0;
}


// ************************************************************************* //
