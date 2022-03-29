/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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

template<class TrackingData>
inline bool Foam::externalPointEdgePoint::update
(
    const point& pt,
    const externalPointEdgePoint& w2,
    const scalar tol,
    TrackingData& td
)
{
    scalar dist2 = magSqr(pt - w2.origin());

    if (!valid(td))
    {
        // current not yet set so use any value
        distSqr_ = dist2;
        origin_ = w2.origin();

        return true;
    }

    scalar diff = distSqr_ - dist2;

    if (diff < 0)
    {
        // already nearer to pt
        return false;
    }

    if ((diff < small) || ((distSqr_ > small) && (diff/distSqr_ < tol)))
    {
        // don't propagate small changes
        return false;
    }
    else
    {
        // update with new values
        distSqr_ = dist2;
        origin_ = w2.origin();

        return true;
    }
}


template<class TrackingData>
inline bool Foam::externalPointEdgePoint::update
(
    const externalPointEdgePoint& w2,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        // current not yet set so use any value
        distSqr_ = w2.distSqr();
        origin_ = w2.origin();

        return true;
    }

    scalar diff = distSqr_ - w2.distSqr();

    if (diff < 0)
    {
        // already nearer to pt
        return false;
    }

    if ((diff < small) || ((distSqr_ > small) && (diff/distSqr_ < tol)))
    {
        // don't propagate small changes
        return false;
    }
    else
    {
        // update with new values
        distSqr_ =  w2.distSqr();
        origin_ = w2.origin();

        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::externalPointEdgePoint::externalPointEdgePoint()
:
    origin_(point::max),
    distSqr_(great)
{}


inline Foam::externalPointEdgePoint::externalPointEdgePoint
(
    const point& origin,
    const scalar distSqr
)
:
    origin_(origin),
    distSqr_(distSqr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point& Foam::externalPointEdgePoint::origin() const
{
    return origin_;
}


inline Foam::scalar Foam::externalPointEdgePoint::distSqr() const
{
    return distSqr_;
}


template<class TrackingData>
inline bool Foam::externalPointEdgePoint::valid(TrackingData& td) const
{
    return origin_ != point::max;
}


// Checks for cyclic points
template<class TrackingData>
inline bool Foam::externalPointEdgePoint::sameGeometry
(
    const externalPointEdgePoint& w2,
    const scalar tol,
    TrackingData& td
) const
{
    scalar diff = Foam::mag(distSqr() - w2.distSqr());

    if (diff < small)
    {
        return true;
    }
    else
    {
        if ((distSqr() > small) && ((diff/distSqr()) < tol))
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
inline void Foam::externalPointEdgePoint::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{
    origin_ = transform.transformPosition(origin_);
}


template<class TrackingData>
inline bool Foam::externalPointEdgePoint::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgeI,
    const externalPointEdgePoint& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(td.points_[pointi], edgeInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::externalPointEdgePoint::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const externalPointEdgePoint& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(td.points_[pointi], newPointInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::externalPointEdgePoint::updatePoint
(
    const externalPointEdgePoint& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(newPointInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::externalPointEdgePoint::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointi,
    const externalPointEdgePoint& pointInfo,
    const scalar tol,
    TrackingData& td
)
{
    const edge& e = mesh.edges()[edgeI];
    return update(e.centre(td.points_), pointInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::externalPointEdgePoint::equal
(
    const externalPointEdgePoint& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::externalPointEdgePoint::operator==
(
    const Foam::externalPointEdgePoint& rhs
)
const
{
    return (origin() == rhs.origin()) && (distSqr() == rhs.distSqr());
}


inline bool Foam::externalPointEdgePoint::operator!=
(
    const Foam::externalPointEdgePoint& rhs
)
const
{
    return !(*this == rhs);
}


// ************************************************************************* //
