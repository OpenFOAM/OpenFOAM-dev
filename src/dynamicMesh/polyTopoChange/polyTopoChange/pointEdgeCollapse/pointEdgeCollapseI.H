/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

// Update this with w2.
template<class TrackingData>
inline bool Foam::pointEdgeCollapse::update
(
    const pointEdgeCollapse& w2,
    const scalar tol,
    TrackingData& td
)
{
    if (!w2.valid(td))
    {
        FatalErrorInFunction
            << "problem." << abort(FatalError);
    }

    if (!valid(td))
    {
        operator=(w2);
        return true;
    }

    if (w2.collapseIndex_ == -1 || collapseIndex_ == -1)
    {
        // Not marked for collapse; only happens on edges.
        return false;
    }

    if (w2.collapsePriority_ < collapsePriority_)
    {
        return false;
    }
    else if (w2.collapsePriority_ > collapsePriority_)
    {
        operator=(w2);
        return true;
    }

    // Get overwritten by w2 if it has a higher priority
    if (w2.collapseIndex_ < collapseIndex_)
    {
        operator=(w2);
        return true;
    }
    else if (w2.collapseIndex_ == collapseIndex_)
    {
        bool identicalPoint = samePoint(w2.collapsePoint_);

        bool nearer = (magSqr(w2.collapsePoint_) < magSqr(collapsePoint_));

        if (nearer)
        {
            operator=(w2);
        }

        if (identicalPoint)
        {
            return false;
        }
        else
        {
            return nearer;
        }
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::pointEdgeCollapse::pointEdgeCollapse()
:
    collapsePoint_(great, great, great),
    collapseIndex_(-2),
    collapsePriority_(-2)
{}


// Construct from origin, distance
inline Foam::pointEdgeCollapse::pointEdgeCollapse
(
    const point& collapsePoint,
    const label collapseIndex,
    const label collapsePriority
)
:
    collapsePoint_(collapsePoint),
    collapseIndex_(collapseIndex),
    collapsePriority_(collapsePriority)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point& Foam::pointEdgeCollapse::collapsePoint() const
{
    return collapsePoint_;
}


inline Foam::label Foam::pointEdgeCollapse::collapseIndex() const
{
    return collapseIndex_;
}


inline Foam::label Foam::pointEdgeCollapse::collapsePriority() const
{
    return collapsePriority_;
}


inline bool Foam::pointEdgeCollapse::samePoint(const point& pt) const
{
    bool isLegal1 = (cmptMin(collapsePoint_) < 0.5*great);
    bool isLegal2 = (cmptMin(pt) < 0.5*great);

    if (isLegal1 && isLegal2)
    {
        return mag(collapsePoint_ - pt) < 1e-9;//small;
    }
    else
    {
        return isLegal1 == isLegal2;
    }
}


template<class TrackingData>
inline bool Foam::pointEdgeCollapse::valid(TrackingData& td) const
{
    return collapseIndex_ != -2;
}


template<class TrackingData>
inline void Foam::pointEdgeCollapse::leaveDomain
(
    const polyPatch& patch,
    const label patchPointi,
    const point& coord,
    TrackingData& td
)
{
    collapsePoint_ -= coord;
}


template<class TrackingData>
inline void Foam::pointEdgeCollapse::transform
(
    const tensor& rotTensor,
    TrackingData& td
)
{
    collapsePoint_ = Foam::transform(rotTensor, collapsePoint_);
}


// Update absolute geometric quantities. Note that distance (dist_)
// is not affected by leaving/entering domain.
template<class TrackingData>
inline void Foam::pointEdgeCollapse::enterDomain
(
    const polyPatch& patch,
    const label patchPointi,
    const point& coord,
    TrackingData& td
)
{
    // back to absolute form
    collapsePoint_ += coord;
}


// Update this with information from connected edge
template<class TrackingData>
inline bool Foam::pointEdgeCollapse::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgeI,
    const pointEdgeCollapse& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(edgeInfo, tol, td);
}


// Update this with new information on same point
template<class TrackingData>
inline bool Foam::pointEdgeCollapse::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const pointEdgeCollapse& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(newPointInfo, tol, td);
}


// Update this with new information on same point. No extra information.
template<class TrackingData>
inline bool Foam::pointEdgeCollapse::updatePoint
(
    const pointEdgeCollapse& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(newPointInfo, tol, td);
}


// Update this with information from connected point
template<class TrackingData>
inline bool Foam::pointEdgeCollapse::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointi,
    const pointEdgeCollapse& pointInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(pointInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::pointEdgeCollapse::equal
(
    const pointEdgeCollapse& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pointEdgeCollapse::operator==
(
    const Foam::pointEdgeCollapse& rhs
) const
{
    return
        collapseIndex_ == rhs.collapseIndex_
     && collapsePriority_ == rhs.collapsePriority_
     && samePoint(rhs.collapsePoint_);
}


inline bool Foam::pointEdgeCollapse::operator!=
(
    const Foam::pointEdgeCollapse& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
