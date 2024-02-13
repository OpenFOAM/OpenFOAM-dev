/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline bool Foam::pistonPointEdgeData::update
(
    const point& pt,
    trackingData& td
)
{
    inBowl_ = (td.axis_& pt) < td.maxZ_;
    const bool propagate = !visited_ && inBowl_;
    visited_ = true;
    return propagate;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::pistonPointEdgeData::pistonPointEdgeData()
:
    visited_(false),
    inBowl_(false)
{}


inline Foam::pistonPointEdgeData::pistonPointEdgeData(const bool inBowl)
:
    visited_(true),
    inBowl_(inBowl)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::pistonPointEdgeData::inBowl() const
{
    return inBowl_;
}


inline bool Foam::pistonPointEdgeData::valid(trackingData& td) const
{
    return visited_;
}


inline void Foam::pistonPointEdgeData::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    trackingData& td
)
{}


inline bool Foam::pistonPointEdgeData::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgei,
    const pistonPointEdgeData& edgeInfo,
    const scalar tol,
    trackingData& td
)
{
    return update(td.points_[pointi], td);
}


inline bool Foam::pistonPointEdgeData::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const pistonPointEdgeData& newPointInfo,
    const scalar tol,
    trackingData& td
)
{
    return update(td.points_[pointi], td);
}


inline bool Foam::pistonPointEdgeData::updateEdge
(
    const polyMesh& mesh,
    const label edgei,
    const label pointi,
    const pistonPointEdgeData& pointInfo,
    const scalar tol,
    trackingData& td
)
{
    visited_ = true;
    return true;
}


inline bool Foam::pistonPointEdgeData::equal
(
    const pistonPointEdgeData& rhs,
    trackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pistonPointEdgeData::operator==
(
    const Foam::pistonPointEdgeData& rhs
)
const
{
    return visited_ == rhs.visited_ && inBowl_ == rhs.inBowl_;
}


inline bool Foam::pistonPointEdgeData::operator!=
(
    const Foam::pistonPointEdgeData& rhs
)
const
{
    return !(*this == rhs);
}


// ************************************************************************* //
