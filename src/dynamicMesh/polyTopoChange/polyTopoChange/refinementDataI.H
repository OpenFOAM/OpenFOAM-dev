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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::refinementData::refinementData()
:
    refinementCount_(-1),
    count_(-1)
{}


// Construct from components
inline Foam::refinementData::refinementData
(
    const label refinementCount,
    const label count
)
:
    refinementCount_(refinementCount),
    count_(count)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::refinementData::valid(TrackingData& td) const
{
    return count_ != -1;
}


// No geometric data so never any problem on cyclics
template<class TrackingData>
inline bool Foam::refinementData::sameGeometry
(
    const polyMesh&,
    const refinementData&,
    const scalar,
    TrackingData& td
) const
{
    return true;
}


// No geometric data.
template<class TrackingData>
inline void Foam::refinementData::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData& td
)
{}


// No geometric data.
template<class TrackingData>
inline void Foam::refinementData::transform
(
    const polyMesh&,
    const tensor& rotTensor,
    TrackingData& td
)
{}


// No geometric data.
template<class TrackingData>
inline void Foam::refinementData::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData& td
)
{}


// Update cell with neighbouring face information
template<class TrackingData>
inline bool Foam::refinementData::updateCell
(
    const polyMesh&,
    const label thisCelli,
    const label neighbourFacei,
    const refinementData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        FatalErrorInFunction
            << abort(FatalError);
        return false;
    }


    // Check if more than 2:1 ratio. This is when I am not refined but neighbour
    // is and neighbour already had higher cell level.
    if
    (
        neighbourInfo.isRefined()
    && !isRefined()
    &&  neighbourInfo.refinementCount() > refinementCount()
    )
    {
        count_ = refinementCount();
        return true;
    }



    // Count from neighbour face by the time it reaches the current cell.
    label transportedFaceCount;

    if (neighbourInfo.isRefined())
    {
        // refined so passes through two cells.
        transportedFaceCount = max(0, neighbourInfo.count()-2);
    }
    else
    {
        // unrefined.
        transportedFaceCount = max(0, neighbourInfo.count()-1);
    }

    if (count_ >= transportedFaceCount)
    {
        return false;
    }
    else
    {
        count_ = transportedFaceCount;

        return true;
    }
}


// Update face with neighbouring cell information
template<class TrackingData>
inline bool Foam::refinementData::updateFace
(
    const polyMesh&,
    const label thisFacei,
    const label neighbourCelli,
    const refinementData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    // From cell to its faces.
    if (!valid(td))
    {
        refinementCount_ = neighbourInfo.refinementCount();
        count_ = neighbourInfo.count();

        return true;
    }

    if (count_ >= neighbourInfo.count())
    {
        return false;
    }
    else
    {
        refinementCount_ = neighbourInfo.refinementCount();
        count_ = neighbourInfo.count();

        return true;
    }
}


// Update face with coupled face information
template<class TrackingData>
inline bool Foam::refinementData::updateFace
(
    const polyMesh&,
    const label thisFacei,
    const refinementData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    // From face to face (e.g. coupled faces)
    if (!valid(td))
    {
        refinementCount_ = neighbourInfo.refinementCount();
        count_ = neighbourInfo.count();

        return true;
    }

    if (count_ >= neighbourInfo.count())
    {
        return false;
    }
    else
    {
        refinementCount_ = neighbourInfo.refinementCount();
        count_ = neighbourInfo.count();

        return true;
    }
}


template<class TrackingData>
inline bool Foam::refinementData::equal
(
    const refinementData& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::refinementData::operator==(const Foam::refinementData& rhs)
 const
{
    return count() == rhs.count() && refinementCount() == rhs.refinementCount();
}


inline bool Foam::refinementData::operator!=(const Foam::refinementData& rhs)
 const
{
    return !(*this == rhs);
}


// ************************************************************************* //
