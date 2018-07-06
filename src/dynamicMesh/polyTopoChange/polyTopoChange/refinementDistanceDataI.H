/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "transform.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Returns the wanted level
inline Foam::label Foam::refinementDistanceData::wantedLevel(const point& pt)
 const
{
    const scalar distSqr = magSqr(pt-origin_);

    // Get the size at the origin level
    scalar levelSize = level0Size_/(1<<originLevel_);

    scalar r = 0;

    for (label level = originLevel_; level >= 0; --level)
    {
        // Current range
        r += levelSize;

        // Check if our distance is within influence sphere
        if (sqr(r) > distSqr)
        {
            return level;
        }

        // Lower level will have double the size
        levelSize *= 2;
    }
    return 0;
}


template<class TrackingData>
inline bool Foam::refinementDistanceData::update
(
    const point& pos,
    const refinementDistanceData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        if (!neighbourInfo.valid(td))
        {
            FatalErrorInFunction
                << "problem" << abort(FatalError);
        }
        operator=(neighbourInfo);
        return true;
    }

    // Determine wanted level at current position.
    label cellLevel = wantedLevel(pos);

    // Determine wanted level coming through the neighbour
    label nbrLevel = neighbourInfo.wantedLevel(pos);

    if (nbrLevel > cellLevel)
    {
        operator=(neighbourInfo);
        return true;
    }
    else if (nbrLevel == cellLevel)
    {
        scalar myDistSqr = magSqr(pos-origin_);
        scalar nbrDistSqr = magSqr(pos - neighbourInfo.origin());
        scalar diff = myDistSqr - nbrDistSqr;

        if (diff < 0)
        {
            // already nearest
            return false;
        }

        if ((diff < small) || ((myDistSqr > small) && (diff/myDistSqr < tol)))
        {
            // don't propagate small changes
            return false;
        }
        else
        {
            // update with new values
            operator=(neighbourInfo);
            return true;
        }
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::refinementDistanceData::refinementDistanceData()
:
    level0Size_(-1)
{}


// Construct from components
inline Foam::refinementDistanceData::refinementDistanceData
(
    const scalar level0Size,
    const point& origin,
    const label originLevel
)
:
    level0Size_(level0Size),
    origin_(origin),
    originLevel_(originLevel)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::refinementDistanceData::valid(TrackingData& td) const
{
    return level0Size_ != -1;
}


// No geometric data so never any problem on cyclics
template<class TrackingData>
inline bool Foam::refinementDistanceData::sameGeometry
(
    const polyMesh&,
    const refinementDistanceData&,
    const scalar,
    TrackingData& td
) const
{
    return true;
}


template<class TrackingData>
inline void Foam::refinementDistanceData::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData& td
)
{
    origin_ -= faceCentre;
}


template<class TrackingData>
inline void Foam::refinementDistanceData::transform
(
    const polyMesh&,
    const tensor& rotTensor,
    TrackingData& td
)
{
    origin_ = Foam::transform(rotTensor, origin_);
}


// Update absolute geometric quantities.
template<class TrackingData>
inline void Foam::refinementDistanceData::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData& td
)
{
    // back to absolute form
    origin_ += faceCentre;
}


// Update cell with neighbouring face information
template<class TrackingData>
inline bool Foam::refinementDistanceData::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label neighbourFacei,
    const refinementDistanceData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    const point& pos = mesh.cellCentres()[thisCelli];

    return update(pos, neighbourInfo, tol, td);
}


// Update face with neighbouring cell information
template<class TrackingData>
inline bool Foam::refinementDistanceData::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const refinementDistanceData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    const point& pos = mesh.faceCentres()[thisFacei];

    return update(pos, neighbourInfo, tol, td);
}


// Update face with coupled face information
template<class TrackingData>
inline bool Foam::refinementDistanceData::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const refinementDistanceData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    const point& pos = mesh.faceCentres()[thisFacei];

    return update(pos, neighbourInfo, tol, td);
}


template<class TrackingData>
inline bool Foam::refinementDistanceData::equal
(
    const refinementDistanceData& rhs,
    TrackingData& td
) const
{
    if (!valid(td))
    {
        if (!rhs.valid(td))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return operator==(rhs);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::refinementDistanceData::operator==
(
    const Foam::refinementDistanceData& rhs
)
 const
{
    return
        level0Size_ == rhs.level0Size_
     && origin_ == rhs.origin_
     && originLevel_ == rhs.originLevel_;
}


inline bool Foam::refinementDistanceData::operator!=
(
    const Foam::refinementDistanceData& rhs
)
 const
{
    return !(*this == rhs);
}


// ************************************************************************* //
