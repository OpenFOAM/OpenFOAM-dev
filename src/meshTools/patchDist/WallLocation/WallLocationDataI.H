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

#include "WallLocationData.H"
#include "transformer.H"

// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

template<class WallLocation, class Type>
template<class TrackingData>
inline bool Foam::WallLocationData<WallLocation, Type>::update
(
    const point& pt,
    const WallLocationData<WallLocation, Type>& w2,
    const scalar tol,
    TrackingData& td
)
{
    const bool result = WallLocation::update(pt, w2, tol, td);

    if (result)
    {
        data_ = w2.data();
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class WallLocation, class Type>
inline Foam::WallLocationData<WallLocation, Type>::WallLocationData()
:
    WallLocation(),
    data_()
{}


template<class WallLocation, class Type>
template<class ... Geometry>
inline Foam::WallLocationData<WallLocation, Type>::WallLocationData
(
    const Type& data,
    const Geometry& ... geometry
)
:
    WallLocation(geometry ...),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class WallLocation, class Type>
inline const Type& Foam::WallLocationData<WallLocation, Type>::data() const
{
    return data_;
}


template<class WallLocation, class Type>
inline Type& Foam::WallLocationData<WallLocation, Type>::data()
{
    return data_;
}


template<class WallLocation, class Type>
template<class TrackingData>
inline const Type& Foam::WallLocationData<WallLocation, Type>::data
(
    TrackingData& td
) const
{
    return data_;
}


template<class WallLocation, class Type>
template<class TrackingData>
inline void Foam::WallLocationData<WallLocation, Type>::transform
(
    const transformer& transform,
    TrackingData& td
)
{
    WallLocation::transform(transform, td);

    data_ = transform.transform(data_);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class WallLocation, class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const WallLocationData<WallLocation, Type>& w
)
{
    return
        os << static_cast<const WallLocation&>(w) << token::SPACE << w.data();
}


template<class WallLocation, class Type>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    WallLocationData<WallLocation, Type>& w
)
{
    return is >> static_cast<WallLocation&>(w) >> w.data_;
}


// ************************************************************************* //
