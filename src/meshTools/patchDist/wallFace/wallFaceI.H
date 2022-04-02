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

#include "wallFace.H"
#include "polyMesh.H"
#include "transformer.H"

// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

template<class Derived>
template<class TrackingData>
inline bool Foam::WallFaceBase<Derived>::update
(
    const point& pt,
    const WallFaceBase<Derived>& w2,
    const scalar tol,
    TrackingData& td
)
{
    const scalar dist2 =
        sqr
        (
            face(identity(w2.points().size()))
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

template<class Derived>
inline Foam::WallFaceBase<Derived>::WallFaceBase()
:
    points_(0),
    distSqr_(-great)
{}


template<class Derived>
inline Foam::WallFaceBase<Derived>::WallFaceBase
(
    const face& f,
    const pointField& points,
    const scalar distSqr
)
:
    points_(f.points(points)),
    distSqr_(distSqr)
{}


template<class Derived>
inline Foam::WallFaceBase<Derived>::WallFaceBase
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

template<class Derived>
inline const Foam::pointField& Foam::WallFaceBase<Derived>::points() const
{
    return points_;
}


template<class Derived>
inline Foam::pointField& Foam::WallFaceBase<Derived>::points()
{
    return points_;
}


template<class Derived>
inline Foam::scalar Foam::WallFaceBase<Derived>::distSqr() const
{
    return distSqr_;
}


template<class Derived>
inline Foam::scalar& Foam::WallFaceBase<Derived>::distSqr()
{
    return distSqr_;
}


template<class Derived>
template<class TrackingData>
inline Foam::scalar Foam::WallFaceBase<Derived>::dist(TrackingData& td) const
{
    return valid(td) ? sqrt(distSqr_) : great;
}


template<class Derived>
template<class TrackingData>
inline bool Foam::WallFaceBase<Derived>::valid(TrackingData& td) const
{
    return distSqr_ > -small;
}


template<class Derived>
template<class TrackingData>
inline bool Foam::WallFaceBase<Derived>::sameGeometry
(
    const polyMesh&,
    const WallFaceBase<Derived>& w2,
    const scalar tol,
    TrackingData& td
) const
{
    scalar diff = mag(distSqr() - w2.distSqr());

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


template<class Derived>
template<class TrackingData>
inline void Foam::WallFaceBase<Derived>::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{
    transform.transformPosition(points_, points_);
}


template<class Derived>
template<class TrackingData>
inline bool Foam::WallFaceBase<Derived>::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label neighbourFacei,
    const WallFaceBase<Derived>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return
        static_cast<Derived&>(*this).update
        (
            mesh.cellCentres()[thisCelli],
            neighbourWallInfo,
            tol,
            td
        );
}


template<class Derived>
template<class TrackingData>
inline bool Foam::WallFaceBase<Derived>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const WallFaceBase<Derived>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return
        static_cast<Derived&>(*this).update
        (
            mesh.faceCentres()[thisFacei],
            neighbourWallInfo,
            tol,
            td
        );
}


template<class Derived>
template<class TrackingData>
inline bool Foam::WallFaceBase<Derived>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const WallFaceBase<Derived>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return
        static_cast<Derived&>(*this).update
        (
            mesh.faceCentres()[thisFacei],
            neighbourWallInfo,
            tol,
            td
        );
}


template<class Derived>
template<class TrackingData>
inline bool Foam::WallFaceBase<Derived>::equal
(
    const WallFaceBase<Derived>& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Derived>
inline bool Foam::WallFaceBase<Derived>::operator==
(
    const Foam::WallFaceBase<Derived>& rhs
) const
{
    return points() == rhs.points();
}


template<class Derived>
inline bool Foam::WallFaceBase<Derived>::operator!=
(
    const Foam::WallFaceBase<Derived>& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Derived>
Foam::Ostream& Foam::operator<<(Ostream& os, const WallFaceBase<Derived>& w)
{
    return os << w.points() << token::SPACE << w.distSqr();
}


template<class Derived>
Foam::Istream& Foam::operator>>(Istream& is, WallFaceBase<Derived>& w)
{
    return is >> w.points() >> w.distSqr();
}


// ************************************************************************* //
