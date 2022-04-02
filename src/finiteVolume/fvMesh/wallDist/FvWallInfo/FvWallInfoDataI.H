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

#include "FvWallInfoData.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class WallInfo, class Type, class Derived>
template<class TrackingData>
inline bool Foam::FvWallInfoDataBase<WallInfo, Type, Derived>::update
(
    const point& pt,
    const FvWallInfoDataBase<WallInfo, Type, Derived>& w2,
    const scalar tol,
    TrackingData& td
)
{
    const bool result =
        FvWallInfoBase<WallInfo, Derived>::update(pt, w2, tol, td);

    if (result)
    {
        data_ = w2.data();
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class WallInfo, class Type, class Derived>
inline Foam::FvWallInfoDataBase<WallInfo, Type, Derived>::FvWallInfoDataBase()
:
    FvWallInfoBase<WallInfo, Derived>(),
    data_()
{}


template<class WallInfo, class Type, class Derived>
template<class ... Geometry>
inline Foam::FvWallInfoDataBase<WallInfo, Type, Derived>::FvWallInfoDataBase
(
    const Type& data,
    const Geometry& ... geometry
)
:
    FvWallInfoBase<WallInfo, Derived>(geometry ...),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class WallInfo, class Type, class Derived>
inline const Type&
Foam::FvWallInfoDataBase<WallInfo, Type, Derived>::data() const
{
    return data_;
}


template<class WallInfo, class Type, class Derived>
inline Type&
Foam::FvWallInfoDataBase<WallInfo, Type, Derived>::data()
{
    return data_;
}


template<class WallInfo, class Type, class Derived>
template<class TrackingData>
inline const Type&
Foam::FvWallInfoDataBase<WallInfo, Type, Derived>::data(TrackingData& td) const
{
    return data_;
}


template<class WallInfo, class Type, class Derived>
template<class TrackingData>
inline void Foam::FvWallInfoDataBase<WallInfo, Type, Derived>::transform
(
    const fvPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{
    FvWallInfoBase<WallInfo, Derived>::transform
    (
        patch,
        patchFacei,
        transform,
        td
    );

    data_ = transform.transform(data_);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class WallInfo, class Type, class Derived>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const FvWallInfoDataBase<WallInfo, Type, Derived>& w
)
{
    return os
        << static_cast<const FvWallInfoBase<WallInfo, Derived>&>(w)
        << token::SPACE
        << w.data();
}


template<class WallInfo, class Type, class Derived>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    FvWallInfoDataBase<WallInfo, Type, Derived>& w
)
{
    return is
        >> static_cast<FvWallInfoBase<WallInfo, Derived>&>(w)
        >> w.data_;
}


// ************************************************************************* //
