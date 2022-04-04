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

#include "smoothDelta.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::LESModels::smoothDelta::deltaData::update
(
    const smoothDelta::deltaData& w2,
    const scalar scale,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td) || (delta_ < vSmall))
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
        // ...

        // Nothing changed
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::LESModels::smoothDelta::deltaData::deltaData()
:
    delta_(-great)
{}


inline Foam::LESModels::smoothDelta::deltaData::deltaData(const scalar delta)
:
    delta_(delta)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::LESModels::smoothDelta::deltaData::valid
(
    TrackingData& td
) const
{
    return delta_ > -small;
}


template<class TrackingData>
inline bool Foam::LESModels::smoothDelta::deltaData::sameGeometry
(
    const fvMesh&,
    const deltaData&,
    const scalar,
    TrackingData& td
) const
{
    return true;
}


template<class TrackingData>
inline void Foam::LESModels::smoothDelta::deltaData::transform
(
    const fvPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


template<class TrackingData>
inline bool Foam::LESModels::smoothDelta::deltaData::updateCell
(
    const fvMesh&,
    const label,
    const labelPair&,
    const deltaData& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    // Take over info from face if more than deltaRatio larger.
    return update(neighbourWallInfo, td, tol, td);
}


template<class TrackingData>
inline bool Foam::LESModels::smoothDelta::deltaData::updateFace
(
    const fvMesh&,
    const labelPair&,
    const label,
    const deltaData& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    // Take over information from cell without any scaling (scale = 1.0)
    return update(neighbourWallInfo, 1.0, tol, td);
}


template<class TrackingData>
inline bool Foam::LESModels::smoothDelta::deltaData::updateFace
(
    const fvMesh&,
    const labelPair&,
    const deltaData& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    // Take over information from coupled face without any scaling (scale = 1.0)
    return update(neighbourWallInfo, 1.0, tol, td);
}


template<class TrackingData>
inline bool Foam::LESModels::smoothDelta::deltaData::equal
(
    const deltaData& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::LESModels::smoothDelta::deltaData::operator==
(
    const deltaData& rhs
) const
{
    return delta_ == rhs.delta();
}


inline bool Foam::LESModels::smoothDelta::deltaData::operator!=
(
    const deltaData& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
