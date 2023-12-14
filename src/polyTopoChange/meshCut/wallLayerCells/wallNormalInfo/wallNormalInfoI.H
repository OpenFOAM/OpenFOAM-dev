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

#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::wallNormalInfo::update
(
    const wallNormalInfo& w2,
    TrackingData& td
)
{
    if (!w2.valid(td))
    {
        FatalErrorInFunction
          << "Problem: w2 is not valid" << abort(FatalError);

        return false;
    }
    else if (valid(td))
    {
        // our already set. Stop any transfer
        return false;
    }
    else
    {
        normal_ = w2.normal();

        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wallNormalInfo::wallNormalInfo()
:
    normal_(point::max)
{}


inline Foam::wallNormalInfo::wallNormalInfo(const vector& normal)
:
    normal_(normal)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::vector& Foam::wallNormalInfo::normal() const
{
    return normal_;
}


template<class TrackingData>
inline bool Foam::wallNormalInfo::valid(TrackingData& td) const
{
    return normal_ != point::max;
}


template<class TrackingData>
inline bool Foam::wallNormalInfo::sameGeometry
(
    const polyMesh&,
    const wallNormalInfo& w2,
    const scalar tol,
    TrackingData& td
) const
{
    return true;
}


template<class TrackingData>
inline void Foam::wallNormalInfo::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


template<class TrackingData>
inline bool Foam::wallNormalInfo::updateCell
(
    const polyMesh&,
    const label thisCelli,
    const label neighbourFacei,
    const wallNormalInfo& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(neighbourWallInfo, td);
}


template<class TrackingData>
inline bool Foam::wallNormalInfo::updateFace
(
    const polyMesh&,
    const label thisFacei,
    const label neighbourCelli,
    const wallNormalInfo& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(neighbourWallInfo, td);
}


template<class TrackingData>
inline bool Foam::wallNormalInfo::updateFace
(
    const polyMesh&,
    const label thisFacei,
    const wallNormalInfo& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update(neighbourWallInfo, td);
}


template<class TrackingData>
inline bool Foam::wallNormalInfo::equal
(
    const wallNormalInfo& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wallNormalInfo::operator==
(
    const Foam::wallNormalInfo& rhs
) const
{
    return normal() == rhs.normal();
}


inline bool Foam::wallNormalInfo::operator!=
(
    const Foam::wallNormalInfo& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
