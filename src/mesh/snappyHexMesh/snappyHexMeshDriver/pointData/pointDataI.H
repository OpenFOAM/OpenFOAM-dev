/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "polyMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::pointData::pointData()
:
    pointEdgePoint(),
    s_(great),
    v_(point::max)
{}


inline Foam::pointData::pointData
(
    const point& origin,
    const scalar distSqr,
    const scalar s,
    const vector& v
)
:
    pointEdgePoint(origin, distSqr),
    s_(s),
    v_(v)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::pointData::s() const
{
    return s_;
}


inline const Foam::vector& Foam::pointData::v() const
{
    return v_;
}


template<class TrackingData>
inline void Foam::pointData::transform
(
    const tensor& rotTensor,
    TrackingData& td
)
{
    pointEdgePoint::transform(rotTensor, td);
    v_ = Foam::transform(rotTensor, v_);
}


// Update this with information from connected edge
template<class TrackingData>
inline bool Foam::pointData::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgeI,
    const pointData& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    if
    (
        pointEdgePoint::updatePoint
        (
            mesh,
            pointi,
            edgeI,
            edgeInfo,
            tol,
            td
        )
    )
    {
        s_ = edgeInfo.s_;
        v_ = edgeInfo.v_;
        return true;
    }
    else
    {
        return false;
    }
}

// Update this with new information on same point
template<class TrackingData>
inline bool Foam::pointData::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const pointData& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if
    (
        pointEdgePoint::updatePoint
        (
            mesh,
            pointi,
            newPointInfo,
            tol,
            td
        )
    )
    {
        s_ = newPointInfo.s_;
        v_ = newPointInfo.v_;
        return true;
    }
    else
    {
        return false;
    }
}


// Update this with new information on same point. No extra information.
template<class TrackingData>
inline bool Foam::pointData::updatePoint
(
    const pointData& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (pointEdgePoint::updatePoint(newPointInfo, tol, td))
    {
        s_ = newPointInfo.s_;
        v_ = newPointInfo.v_;
        return true;
    }
    else
    {
        return false;
    }
}


// Update this with information from connected point
template<class TrackingData>
inline bool Foam::pointData::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointi,
    const pointData& pointInfo,
    const scalar tol,
    TrackingData& td

)
{
    if
    (
        pointEdgePoint::updateEdge
        (
            mesh,
            edgeI,
            pointi,
            pointInfo,
            tol,
            td
        )
    )
    {
        s_ = pointInfo.s_;
        v_ = pointInfo.v_;
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pointData::operator==(const Foam::pointData& rhs)
const
{
    return
        pointEdgePoint::operator==(rhs)
     && (s() == rhs.s())
     && (v() == rhs.v());
}


inline bool Foam::pointData::operator!=(const Foam::pointData& rhs)
const
{
    return !(*this == rhs);
}


// ************************************************************************* //
