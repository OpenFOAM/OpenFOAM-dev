/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2026 OpenFOAM Foundation
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

#include "pointTopoDistanceData.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::pointTopoDistanceData::pointTopoDistanceData()
:
    data_(-1),
    distance_(-1)
{}


inline Foam::pointTopoDistanceData::pointTopoDistanceData
(
    const label data,
    const label distance
)
:
    data_(data),
    distance_(distance)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::pointTopoDistanceData::valid(TrackingData& td) const
{
    return distance_ != -1;
}


template<class TrackingData>
inline bool Foam::pointTopoDistanceData::sameGeometry
(
    const pointTopoDistanceData&,
    const scalar tol,
    TrackingData& td
) const
{
    return true;
}


template<class TrackingData>
inline void Foam::pointTopoDistanceData::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


template<class TrackingData>
inline bool Foam::pointTopoDistanceData::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgeI,
    const pointTopoDistanceData& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        data_ = edgeInfo.data_;
        distance_ = edgeInfo.distance_ + 1;
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::pointTopoDistanceData::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const pointTopoDistanceData& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        operator=(newPointInfo);
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::pointTopoDistanceData::updatePoint
(
    const pointTopoDistanceData& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        operator=(newPointInfo);
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::pointTopoDistanceData::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointi,
    const pointTopoDistanceData& pointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        operator=(pointInfo);
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::pointTopoDistanceData::equal
(
    const pointTopoDistanceData& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pointTopoDistanceData::operator==
(
    const Foam::pointTopoDistanceData& rhs
) const
{
    return data() == rhs.data() && distance() == rhs.distance();
}


inline bool Foam::pointTopoDistanceData::operator!=
(
    const Foam::pointTopoDistanceData& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::pointTopoDistanceData& wDist
)
{
    return os << wDist.data_ << token::SPACE << wDist.distance_;
}


inline Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::pointTopoDistanceData& wDist
)
{
    return is >> wDist.data_ >> wDist.distance_;
}


// ************************************************************************* //
