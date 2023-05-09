/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "PointEdgeLayerInfoData.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::PointEdgeLayerInfoData<Type>::PointEdgeLayerInfoData()
:
    pointEdgeLayerInfo(),
    data_()
{}


template<class Type>
inline Foam::PointEdgeLayerInfoData<Type>::PointEdgeLayerInfoData
(
    const label pointLayer,
    const Type& data
)
:
    pointEdgeLayerInfo(pointLayer),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline const Type& Foam::PointEdgeLayerInfoData<Type>::data() const
{
    return data_;
}


template<class Type>
template<class TrackingData>
inline void Foam::PointEdgeLayerInfoData<Type>::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{
    if (transform.transformsPosition())
    {
        FatalErrorInFunction
            << "Cannot have a transform within a set of mesh layers"
            << exit(FatalError);
    }
}


template<class Type>
template<class TrackingData>
inline bool Foam::PointEdgeLayerInfoData<Type>::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgei,
    const PointEdgeLayerInfoData<Type>& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (pointEdgeLayerInfo::updatePoint(mesh, pointi, edgei, edgeInfo, tol, td))
    {
        data_ = edgeInfo.data_;
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
template<class TrackingData>
inline bool Foam::PointEdgeLayerInfoData<Type>::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const PointEdgeLayerInfoData<Type>& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (pointEdgeLayerInfo::updatePoint(mesh, pointi, newPointInfo, tol, td))
    {
        data_ = newPointInfo.data_;
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
template<class TrackingData>
inline bool Foam::PointEdgeLayerInfoData<Type>::updatePoint
(
    const PointEdgeLayerInfoData<Type>& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (pointEdgeLayerInfo::updatePoint(newPointInfo, tol, td))
    {
        data_ = newPointInfo.data_;
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
template<class TrackingData>
inline bool Foam::PointEdgeLayerInfoData<Type>::updateEdge
(
    const polyMesh& mesh,
    const label edgei,
    const label pointi,
    const PointEdgeLayerInfoData<Type>& pointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (pointEdgeLayerInfo::updateEdge(mesh, edgei, pointi, pointInfo, tol, td))
    {
        data_ = pointInfo.data_;
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
template<class TrackingData>
inline bool Foam::PointEdgeLayerInfoData<Type>::equal
(
    const PointEdgeLayerInfoData<Type>& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::PointEdgeLayerInfoData<Type>::operator==
(
    const Foam::PointEdgeLayerInfoData<Type>& rhs
) const
{
    return pointEdgeLayerInfo::operator==(rhs) && data_ == rhs.data_;
}


template<class Type>
inline bool Foam::PointEdgeLayerInfoData<Type>::operator!=
(
    const Foam::PointEdgeLayerInfoData<Type>& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PointEdgeLayerInfoData<Type>& l
)
{
    return
        os
     << static_cast<const pointEdgeLayerInfo&>(l)
     << token::SPACE
     << l.data_;
}


template<class Type>
Foam::Istream& Foam::operator>>(Istream& is, PointEdgeLayerInfoData<Type>& l)
{
    return is >> static_cast<pointEdgeLayerInfo&>(l) >> l.data_;
}


// ************************************************************************* //
