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

#include "PatchEdgeFacePointData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::PatchEdgeFacePointData<Type>::PatchEdgeFacePointData()
:
    patchEdgeFacePoint(),
    data_()
{}


template<class Type>
inline Foam::PatchEdgeFacePointData<Type>::PatchEdgeFacePointData
(
    const Type& data,
    const point& origin,
    const scalar distSqr
)
:
    patchEdgeFacePoint(origin, distSqr),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline const Type& Foam::PatchEdgeFacePointData<Type>::data() const
{
    return data_;
}


template<class Type>
inline Type& Foam::PatchEdgeFacePointData<Type>::data()
{
    return data_;
}


template<class Type>
template<class TrackingData>
inline const Type& Foam::PatchEdgeFacePointData<Type>::data
(
    TrackingData& td
) const
{
    return data_;
}


template<class Type>
template<class TrackingData>
inline void Foam::PatchEdgeFacePointData<Type>::transform
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const tensor& rotTensor,
    const scalar tol,
    TrackingData& td
)
{
    patchEdgeFacePoint::transform(mesh, patch, rotTensor, tol, td);

    data_ = Foam::transform(rotTensor, data_);
}


template<class Type>
template<class TrackingData>
inline bool Foam::PatchEdgeFacePointData<Type>::updateEdge
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const label edgei,
    const label facei,
    const PatchEdgeFacePointData<Type>& faceInfo,
    const scalar tol,
    TrackingData& td
)
{
    const bool result =
        patchEdgeFacePoint::updateEdge
        (
            mesh,
            patch,
            edgei,
            facei,
            faceInfo,
            tol,
            td
        );

    if (result)
    {
        data_ = faceInfo.data_;
    }

    return result;
}


template<class Type>
template<class TrackingData>
inline bool Foam::PatchEdgeFacePointData<Type>::updateEdge
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const PatchEdgeFacePointData<Type>& edgeInfo,
    const bool sameOrientation,
    const scalar tol,
    TrackingData& td
)
{
    const bool result =
        patchEdgeFacePoint::updateEdge
        (
            mesh,
            patch,
            edgeInfo,
            sameOrientation,
            tol,
            td
        );

    if (result)
    {
        data_ = edgeInfo.data_;
    }

    return result;
}


template<class Type>
template<class TrackingData>
inline bool Foam::PatchEdgeFacePointData<Type>::updateFace
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const label facei,
    const label edgei,
    const PatchEdgeFacePointData<Type>& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    const bool result =
        patchEdgeFacePoint::updateFace
        (
            mesh,
            patch,
            facei,
            edgei,
            edgeInfo,
            tol,
            td
        );

    if (result)
    {
        data_ = edgeInfo.data_;
    }

    return result;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PatchEdgeFacePointData<Type>& wDist
)
{
    return
        os
     << static_cast<const patchEdgeFacePoint&>(wDist)
     << token::SPACE
     << wDist.data_;
}


template<class Type>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    PatchEdgeFacePointData<Type>& wDist
)
{
    return is >> static_cast<patchEdgeFacePoint&>(wDist) >> wDist.data_;
}


// ************************************************************************* //
