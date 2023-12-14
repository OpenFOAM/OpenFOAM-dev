/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

inline bool Foam::pointEdgeDist::update
(
    const point& pt,
    const pointEdgeDist& w2,
    const scalar tol,
    data& td
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::pointEdgeDist::pointEdgeDist()
:
    origin_(point::max),
    distSqr_(great)
{}


inline Foam::pointEdgeDist::pointEdgeDist
(
    const point& origin,
    const scalar distSqr
)
:
    origin_(origin),
    distSqr_(distSqr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point& Foam::pointEdgeDist::origin() const
{
    return origin_;
}


inline Foam::scalar Foam::pointEdgeDist::distSqr() const
{
    return distSqr_;
}


inline bool Foam::pointEdgeDist::valid(data& td) const
{
    return origin_ != point::max;
}


inline void Foam::pointEdgeDist::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    data& td
)
{
    origin_ = transform.transformPosition(origin_);
}


inline bool Foam::pointEdgeDist::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgei,
    const pointEdgeDist& edgeInfo,
    const scalar tol,
    data& td
)
{
    return update(td.points[pointi], edgeInfo, tol, td);
}


inline bool Foam::pointEdgeDist::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const pointEdgeDist& newPointInfo,
    const scalar tol,
    data& td
)
{
    return update(td.points[pointi], newPointInfo, tol, td);
}


inline bool Foam::pointEdgeDist::updateEdge
(
    const polyMesh& mesh,
    const label edgei,
    const label pointi,
    const pointEdgeDist& pointInfo,
    const scalar tol,
    data& td
)
{
    const edge& e = mesh.edges()[edgei];
    return update(e.centre(td.points), pointInfo, tol, td);
}


inline bool Foam::pointEdgeDist::equal
(
    const pointEdgeDist& rhs,
    data& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pointEdgeDist::operator==
(
    const Foam::pointEdgeDist& rhs
)
const
{
    return (origin() == rhs.origin()) && (distSqr() == rhs.distSqr());
}


inline bool Foam::pointEdgeDist::operator!=
(
    const Foam::pointEdgeDist& rhs
)
const
{
    return !(*this == rhs);
}


// ************************************************************************* //
