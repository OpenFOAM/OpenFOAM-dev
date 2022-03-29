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

#include "polyMesh.H"
#include "transformer.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2.
template<class TrackingData>
inline bool Foam::pointEdgeStructuredWalk::update
(
    const pointEdgeStructuredWalk& w2,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        // current not yet set. Walked from w2 to here (=point0)
        dist_ = w2.dist_ + mag(point0_-w2.previousPoint_);
        previousPoint_ = point0_;
        data_ = w2.data_;

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::pointEdgeStructuredWalk::pointEdgeStructuredWalk()
:
    point0_(vector::max),
    previousPoint_(vector::max),
    dist_(0),
    data_(Zero)
{}


// Construct from origin, distance
inline Foam::pointEdgeStructuredWalk::pointEdgeStructuredWalk
(
    const point& point0,
    const point& previousPoint,
    const scalar dist,
    const vector& data
)
:
    point0_(point0),
    previousPoint_(previousPoint),
    dist_(dist),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::pointEdgeStructuredWalk::inZone() const
{
    return point0_ != vector::max;
}


//inline const Foam::point& Foam::pointEdgeStructuredWalk::previousPoint() const
//{
//    return previousPoint_;
//}


inline Foam::scalar Foam::pointEdgeStructuredWalk::dist() const
{
    return dist_;
}


inline const Foam::vector& Foam::pointEdgeStructuredWalk::data() const
{
    return data_;
}


template<class TrackingData>
inline bool Foam::pointEdgeStructuredWalk::valid(TrackingData& td) const
{
    return previousPoint_ != vector::max;
}


// Checks for cyclic points
template<class TrackingData>
inline bool Foam::pointEdgeStructuredWalk::sameGeometry
(
    const pointEdgeStructuredWalk& w2,
    const scalar tol,
    TrackingData& td
) const
{
    scalar diff = Foam::mag(dist() - w2.dist());

    if (diff < small)
    {
        return true;
    }
    else
    {
        if ((dist() > small) && ((diff/dist()) < tol))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


template<class TrackingData>
inline void Foam::pointEdgeStructuredWalk::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{
    // Note that dist_ is not affected by crossing an interface
    previousPoint_ = transform.transformPosition(previousPoint_);
}


// Update this with information from connected edge
template<class TrackingData>
inline bool Foam::pointEdgeStructuredWalk::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgeI,
    const pointEdgeStructuredWalk& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (inZone())
    {
        return update(edgeInfo, tol, td);
    }
    else
    {
        return false;
    }
}


// Update this with new information on same point
template<class TrackingData>
inline bool Foam::pointEdgeStructuredWalk::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const pointEdgeStructuredWalk& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (inZone())
    {
        return update(newPointInfo, tol, td);
    }
    else
    {
        return false;
    }
}


// Update this with new information on same point. No extra information.
template<class TrackingData>
inline bool Foam::pointEdgeStructuredWalk::updatePoint
(
    const pointEdgeStructuredWalk& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(newPointInfo, tol, td);
}


// Update this with information from connected point
template<class TrackingData>
inline bool Foam::pointEdgeStructuredWalk::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointi,
    const pointEdgeStructuredWalk& pointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (inZone())
    {
        return update(pointInfo, tol, td);
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::pointEdgeStructuredWalk::equal
(
    const pointEdgeStructuredWalk& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pointEdgeStructuredWalk::operator==
(
    const Foam::pointEdgeStructuredWalk& rhs
) const
{
    return previousPoint_ == rhs.previousPoint_;
}


inline bool Foam::pointEdgeStructuredWalk::operator!=
(
    const Foam::pointEdgeStructuredWalk& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
