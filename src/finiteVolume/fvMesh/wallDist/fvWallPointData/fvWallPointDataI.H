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

#include "fvWallPointData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class TrackingData>
inline bool fvWallPointData<Type>::update
(
    const point& pt,
    const fvWallPointData<Type>& w2,
    const scalar tol,
    TrackingData& td
)
{
    scalar dist2 = magSqr(pt - w2.origin());

    if (valid(td))
    {
        scalar diff = distSqr() - dist2;

        if (diff < 0)
        {
            // already nearer to pt
            return false;
        }

        if ((diff < small) || ((distSqr() > small) && (diff/distSqr() < tol)))
        {
            // don't propagate small changes
            return false;
        }
    }

    // Either *this is not yet valid or w2 is closer
    {
        // current not yet set so use any value
        distSqr() = dist2;
        origin() = w2.origin();
        data_ = w2.data();

        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline fvWallPointData<Type>::fvWallPointData()
:
    fvWallPoint(),
    data_()
{}


template<class Type>
inline fvWallPointData<Type>::fvWallPointData
(
    const point& origin,
    const Type& data,
    const scalar distSqr
)
:
    fvWallPoint(origin, distSqr),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline const Type& fvWallPointData<Type>::data() const
{
    return data_;
}


template<class Type>
inline Type& fvWallPointData<Type>::data()
{
    return data_;
}


template<class Type>
template<class TrackingData>
inline const Type& fvWallPointData<Type>::data(TrackingData& td) const
{
    return data_;
}


template<class Type>
template<class TrackingData>
inline void Foam::fvWallPointData<Type>::transform
(
    const fvPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{
    fvWallPoint::transform(patch, patchFacei, transform, td);
    data_ = transform.transform(data_);
}


template<class Type>
template<class TrackingData>
inline bool fvWallPointData<Type>::updateCell
(
    const fvMesh& mesh,
    const label thisCelli,
    const labelPair& neighbourPatchAndFacei,
    const fvWallPointData<Type>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update
    (
        mesh.C()[thisCelli],
        neighbourWallInfo,
        tol,
        td
    );
}


template<class Type>
template<class TrackingData>
inline bool fvWallPointData<Type>::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const label neighbourCelli,
    const fvWallPointData<Type>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label thisPatchi = thisPatchAndFacei.first();
    const label thisFacei = thisPatchAndFacei.second();

    return update
    (
        thisPatchi == -1
      ? mesh.Cf()[thisFacei]
      : mesh.Cf().boundaryField()[thisPatchi][thisFacei],
        neighbourWallInfo,
        tol,
        td
    );
}


template<class Type>
template<class TrackingData>
inline bool fvWallPointData<Type>::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const fvWallPointData<Type>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label thisPatchi = thisPatchAndFacei.first();
    const label thisFacei = thisPatchAndFacei.second();

    return update
    (
        thisPatchi == -1
      ? mesh.Cf()[thisFacei]
      : mesh.Cf().boundaryField()[thisPatchi][thisFacei],
        neighbourWallInfo,
        tol,
        td
    );
}


// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
