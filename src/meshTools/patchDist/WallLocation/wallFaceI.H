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

#include "wallFace.H"
#include "polyMesh.H"
#include "transformer.H"

// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::wallFace::update
(
    const point& pt,
    const wallFace& w2,
    const scalar tol,
    TrackingData& td
)
{
    const scalar dist2 =
        sqr
        (
            face(identityMap(w2.points().size()))
           .nearestPoint(pt, w2.points())
           .distance()
        );

    if (valid(td))
    {
        scalar diff = distSqr() - dist2;

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
    points() = w2.points();

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wallFace::wallFace()
:
    points_(0),
    distSqr_(-great)
{}


inline Foam::wallFace::wallFace
(
    const face& f,
    const pointField& points,
    const scalar distSqr
)
:
    points_(f.points(points)),
    distSqr_(distSqr)
{}


inline Foam::wallFace::wallFace
(
    const face& f,
    const pointField& points,
    const point& centre,
    const scalar distSqr
)
:
    points_(f.points(points)),
    distSqr_(distSqr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::pointField& Foam::wallFace::points() const
{
    return points_;
}


inline Foam::pointField& Foam::wallFace::points()
{
    return points_;
}


inline Foam::scalar Foam::wallFace::distSqr() const
{
    return distSqr_;
}


inline Foam::scalar& Foam::wallFace::distSqr()
{
    return distSqr_;
}


template<class TrackingData>
inline Foam::scalar Foam::wallFace::dist(TrackingData& td) const
{
    return valid(td) ? sqrt(distSqr_) : great;
}


template<class TrackingData>
inline bool Foam::wallFace::valid(TrackingData& td) const
{
    return distSqr_ > -small;
}


template<class TrackingData>
inline bool Foam::wallFace::sameGeometry
(
    const wallFace& w2,
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
inline void Foam::wallFace::transform
(
    const transformer& transform,
    TrackingData& td
)
{
    // Note that distSqr_ is not affected by crossing an interface
    transform.transformPosition(points_, points_);
}


template<class TrackingData>
inline bool Foam::wallFace::equal
(
    const wallFace& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wallFace::operator==
(
    const Foam::wallFace& rhs
) const
{
    return points() == rhs.points();
}


inline bool Foam::wallFace::operator!=
(
    const Foam::wallFace& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const wallFace& w)
{
    return os << w.points() << token::SPACE << w.distSqr();
}


inline Foam::Istream& Foam::operator>>(Istream& is, wallFace& w)
{
    return is >> w.points() >> w.distSqr();
}


// ************************************************************************* //
