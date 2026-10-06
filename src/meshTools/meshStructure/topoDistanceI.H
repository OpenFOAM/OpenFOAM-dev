/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "topoDistance.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::topoDistance::topoDistance()
:
    distance_(-1)
{}


inline Foam::topoDistance::topoDistance(const label distance)
:
    distance_(distance)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::topoDistance::valid(TrackingData& td) const
{
    return distance_ != -1;
}


template<class TrackingData>
inline bool Foam::topoDistance::sameGeometry
(
    const polyMesh&,
    const topoDistance&,
    const scalar,
    TrackingData&
) const
{
    return true;
}


template<class Patch, class TrackingData>
inline void Foam::topoDistance::transform
(
    const Patch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


template<class Mesh, class NeighbourFacei, class TrackingData>
inline bool Foam::topoDistance::updateCell
(
    const Mesh&,
    const label thisCelli,
    const NeighbourFacei& neighbourFacei,
    const topoDistance& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    if (distance_ == -1)
    {
        operator=(neighbourInfo);
        return true;
    }
    else
    {
        return false;
    }
}


template<class Mesh, class ThisFacei, class TrackingData>
inline bool Foam::topoDistance::updateFace
(
    const Mesh&,
    const ThisFacei& thisFacei,
    const label neighbourCelli,
    const topoDistance& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    if (distance_ == -1)
    {
        distance_ = neighbourInfo.distance_ + 1;
        return true;
    }
    else
    {
        return false;
    }
}


template<class Mesh, class ThisFacei, class TrackingData>
inline bool Foam::topoDistance::updateFace
(
    const Mesh&,
    const ThisFacei& thisFacei,
    const topoDistance& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    if (distance_ == -1)
    {
        operator=(neighbourInfo);
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::topoDistance::equal
(
    const topoDistance& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::topoDistance::operator==
(
    const Foam::topoDistance& rhs
) const
{
    return distance() == rhs.distance();
}


inline bool Foam::topoDistance::operator!=
(
    const Foam::topoDistance& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::topoDistance& wDist
)
{
    return os << wDist.distance_;
}


inline Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::topoDistance& wDist
)
{
    return is >> wDist.distance_;
}


// ************************************************************************* //
