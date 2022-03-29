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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::topoDistanceData::topoDistanceData()
:
    data_(-1),
    distance_(-1)
{}


// Construct from components
inline Foam::topoDistanceData::topoDistanceData
(
    const label data,
    const label distance
)
:
    data_(data),
    distance_(distance)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::topoDistanceData::valid(TrackingData& td) const
{
    return distance_ != -1;
}


// No geometric data so never any problem on cyclics
template<class TrackingData>
inline bool Foam::topoDistanceData::sameGeometry
(
    const polyMesh&,
    const topoDistanceData&,
    const scalar,
    TrackingData&
) const
{
    return true;
}


// No geometric data.
template<class TrackingData>
inline void Foam::topoDistanceData::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


// Update cell with neighbouring face information
template<class TrackingData>
inline bool Foam::topoDistanceData::updateCell
(
    const polyMesh&,
    const label thisCelli,
    const label neighbourFacei,
    const topoDistanceData& neighbourInfo,
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


// Update face with neighbouring cell information
template<class TrackingData>
inline bool Foam::topoDistanceData::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const topoDistanceData& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    // From cell to its faces.

    if (distance_ == -1)
    {
        data_ = neighbourInfo.data_;
        distance_ = neighbourInfo.distance_ + 1;
        return true;
    }
    else
    {
        return false;
    }
}


// Update face with coupled face information
template<class TrackingData>
inline bool Foam::topoDistanceData::updateFace
(
    const polyMesh&,
    const label thisFacei,
    const topoDistanceData& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    // From face to face (e.g. coupled faces)
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
inline bool Foam::topoDistanceData::equal
(
    const topoDistanceData& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::topoDistanceData::operator==
(
    const Foam::topoDistanceData& rhs
) const
{
    return data() == rhs.data() && distance() == rhs.distance();
}


inline bool Foam::topoDistanceData::operator!=
(
    const Foam::topoDistanceData& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
