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

#include "wallPoint.H"
#include "transformer.H"
#include "SubField.H"

// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::wallPoint::update
(
    const point& pt,
    const wallPoint& w2,
    const scalar tol,
    TrackingData& td
)
{
    const scalar dist2 = magSqr(pt - w2.origin());

    if (valid(td))
    {
        const scalar diff = distSqr() - dist2;

        if (diff < 0)
        {
            // already nearer to pt
            return false;
        }

        if ((diff < small) || ((distSqr() > small) && (diff/distSqr() < tol)))
        {
            // don't propagate small changes
            return false;
        }
    }

    // Either *this is not yet valid or w2 is closer
    distSqr() = dist2;
    origin() = w2.origin();

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wallPoint::wallPoint()
:
    origin_(point::max),
    distSqr_(-great)
{}


inline Foam::wallPoint::wallPoint
(
    const point& origin,
    const scalar distSqr
)
:
    origin_(origin),
    distSqr_(distSqr)
{}


inline Foam::wallPoint::wallPoint
(
    const face& f,
    const pointField& points,
    const point& centre,
    const scalar distSqr
)
:
    origin_(centre),
    distSqr_(distSqr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point& Foam::wallPoint::origin() const
{
    return origin_;
}


inline Foam::point& Foam::wallPoint::origin()
{
    return origin_;
}


inline Foam::scalar Foam::wallPoint::distSqr() const
{
    return distSqr_;
}


inline Foam::scalar& Foam::wallPoint::distSqr()
{
    return distSqr_;
}


template<class TrackingData>
inline Foam::scalar Foam::wallPoint::dist(TrackingData& td) const
{
    return valid(td) ? sqrt(distSqr_) : great;
}


template<class TrackingData>
inline bool Foam::wallPoint::valid(TrackingData& td) const
{
    return distSqr_ > -small;
}


template<class TrackingData>
inline bool Foam::wallPoint::sameGeometry
(
    const wallPoint& w2,
    const scalar tol,
    TrackingData& td
) const
{
    const scalar diff = mag(distSqr() - w2.distSqr());

    return
        diff < small
     || ((distSqr() > small) && (diff/distSqr() < tol));
}


template<class TrackingData>
inline void Foam::wallPoint::transform
(
    const transformer& transform,
    TrackingData& td
)
{
    // Note that distSqr_ is not affected by crossing an interface
    origin_ = transform.transformPosition(origin_);
}


template<class TrackingData>
inline bool Foam::wallPoint::equal
(
    const wallPoint& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wallPoint::operator==
(
    const Foam::wallPoint& rhs
) const
{
    return origin() == rhs.origin();
}


inline bool Foam::wallPoint::operator!=
(
    const Foam::wallPoint& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const wallPoint& w)
{
    return os << w.origin() << token::SPACE << w.distSqr();
}


inline Foam::Istream& Foam::operator>>(Istream& is, wallPoint& w)
{
    return is >> w.origin() >> w.distSqr();
}


// ************************************************************************* //
