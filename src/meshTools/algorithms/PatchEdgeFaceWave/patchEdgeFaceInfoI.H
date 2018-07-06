/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::patchEdgeFaceInfo::update
(
    const point& pt,
    const patchEdgeFaceInfo& w2,
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


// Update this with w2 (information on same edge)
template<class TrackingData>
inline bool Foam::patchEdgeFaceInfo::update
(
    const patchEdgeFaceInfo& w2,
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

// Null constructor
inline Foam::patchEdgeFaceInfo::patchEdgeFaceInfo()
:
    origin_(point::max),
    distSqr_(sqr(great))
{}


// Construct from origin, distance
inline Foam::patchEdgeFaceInfo::patchEdgeFaceInfo
(
    const point& origin,
    const scalar distSqr
)
:
    origin_(origin),
    distSqr_(distSqr)
{}


// Construct as copy
inline Foam::patchEdgeFaceInfo::patchEdgeFaceInfo(const patchEdgeFaceInfo& wpt)
:
    origin_(wpt.origin()),
    distSqr_(wpt.distSqr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point& Foam::patchEdgeFaceInfo::origin() const
{
    return origin_;
}


inline Foam::scalar Foam::patchEdgeFaceInfo::distSqr() const
{
    return distSqr_;
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceInfo::valid(TrackingData& td) const
{
    return origin_ != point::max;
}


template<class TrackingData>
inline void Foam::patchEdgeFaceInfo::transform
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const tensor& rotTensor,
    const scalar tol,
    TrackingData& td
)
{
    origin_ = Foam::transform(rotTensor, origin_);
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceInfo::updateEdge
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const label edgeI,
    const label facei,
    const patchEdgeFaceInfo& faceInfo,
    const scalar tol,
    TrackingData& td
)
{
    const edge& e = patch.edges()[edgeI];
    point eMid =
        0.5
      * (
            patch.points()[patch.meshPoints()[e[0]]]
          + patch.points()[patch.meshPoints()[e[1]]]
        );
    return update(eMid, faceInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceInfo::updateEdge
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const patchEdgeFaceInfo& edgeInfo,
    const bool sameOrientation,
    const scalar tol,
    TrackingData& td
)
{
    return update(edgeInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceInfo::updateFace
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const label facei,
    const label edgeI,
    const patchEdgeFaceInfo& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(patch.faceCentres()[facei], edgeInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceInfo::equal
(
    const patchEdgeFaceInfo& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::patchEdgeFaceInfo::operator==
(
    const Foam::patchEdgeFaceInfo& rhs
) const
{
    return origin() == rhs.origin();
}


inline bool Foam::patchEdgeFaceInfo::operator!=
(
    const Foam::patchEdgeFaceInfo& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
