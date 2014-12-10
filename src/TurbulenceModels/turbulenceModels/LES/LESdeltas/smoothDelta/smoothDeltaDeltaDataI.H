/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 if applicable
template<class TrackingData>
inline bool smoothDelta::deltaData::update
(
    const smoothDelta::deltaData& w2,
    const scalar scale,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td) || (delta_ < VSMALL))
    {
        // My delta not set. Take over neighbour.
        delta_ = w2.delta()/scale;

        // Something changed. Let caller know.
        return true;
    }
    else if (w2.delta() > (1 + tol)*scale*delta_)
    {
        // Neighbour is too big for me. Up my delta.
        delta_ = w2.delta()/scale;

        // Something changed. Let caller know.
        return true;
    }
    else
    {
        // Neighbour is not too big for me or change is too small
        // Nothing changed.
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline smoothDelta::deltaData::deltaData()
:
    delta_(-GREAT)
{}


// Construct from components
inline smoothDelta::deltaData::deltaData(const scalar delta)
:
    delta_(delta)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool smoothDelta::deltaData::valid(TrackingData& td) const
{
    return delta_ > -SMALL;
}


// Checks for cyclic faces
template<class TrackingData>
inline bool smoothDelta::deltaData::sameGeometry
(
    const polyMesh&,
    const deltaData&,
    const scalar,
    TrackingData& td
) const
{
    return true;
}


template<class TrackingData>
inline void smoothDelta::deltaData::leaveDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point&,
    TrackingData& td
)
{}


template<class TrackingData>
inline void smoothDelta::deltaData::transform
(
    const polyMesh&,
    const tensor&,
    TrackingData& td
)
{}


// Update absolute geometric quantities.
template<class TrackingData>
inline void smoothDelta::deltaData::enterDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point&,
    TrackingData& td
)
{}


// Update this (cellI) with face information.
template<class TrackingData>
inline bool smoothDelta::deltaData::updateCell
(
    const polyMesh&,
    const label,
    const label,
    const deltaData& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    // Take over info from face if more than deltaRatio larger.
    return update(neighbourWallInfo, td, tol, td);
}


// Update this (face) with cell information.
template<class TrackingData>
inline bool smoothDelta::deltaData::updateFace
(
    const polyMesh&,
    const label,
    const label,
    const deltaData& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    // Take over information from cell without any scaling (scale = 1.0)
    return update(neighbourWallInfo, 1.0, tol, td);
}


// Update this (face) with coupled face information.
template<class TrackingData>
inline bool smoothDelta::deltaData::updateFace
(
    const polyMesh&,
    const label,
    const deltaData& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    // Take over information from coupled face without any scaling (scale = 1.0)
    return update(neighbourWallInfo, 1.0, tol, td);
}


template<class TrackingData>
inline bool smoothDelta::deltaData::equal
(
    const deltaData& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool smoothDelta::deltaData::operator==
(
    const deltaData& rhs
) const
{
    return delta_ == rhs.delta();
}


inline bool smoothDelta::deltaData::operator!=
(
    const deltaData& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
