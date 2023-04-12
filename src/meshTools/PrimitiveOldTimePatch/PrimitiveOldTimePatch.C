/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "PrimitiveOldTimePatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::PrimitiveOldTimePatch<FaceList, PointField>::PrimitiveOldTimePatch
(
    const FaceList& faces,
    const Field<PointType>& points,
    const Field<PointType>& points0
)
:
    PrimitivePatch<FaceList, PointField>(faces, points),
    points0Ptr_(isRef<PointField> ? nullptr : new Field<PointType>(points0)),
    points0_(isRef<PointField> ? points0 : points0Ptr_()),
    patch0Ptr_(new patch0Type(faces, points0_)),
    localPoints0Ptr_(nullptr)
{}


template<class FaceList, class PointField>
Foam::PrimitiveOldTimePatch<FaceList, PointField>::PrimitiveOldTimePatch
(
    const PrimitivePatch<FaceList, PointField>& patch,
    const Field<PointType>& points0
)
:
    PrimitivePatch<FaceList, PointField>(patch),
    points0_(points0),
    patch0Ptr_(new patch0Type(patch, points0_)),
    localPoints0Ptr_(nullptr)
{}


template<class FaceList, class PointField>
Foam::PrimitiveOldTimePatch<FaceList, PointField>::PrimitiveOldTimePatch
(
    const FaceList& faces,
    const Field<PointType>& points
)
:
    PrimitivePatch<FaceList, PointField>(faces, points),
    points0Ptr_(nullptr),
    points0_(NullObjectRef<Field<PointType>>()),
    patch0Ptr_(nullptr),
    localPoints0Ptr_(nullptr)
{}


template<class FaceList, class PointField>
Foam::PrimitiveOldTimePatch<FaceList, PointField>::PrimitiveOldTimePatch
(
    const PrimitivePatch<FaceList, PointField>& patch
)
:
    PrimitivePatch<FaceList, PointField>(patch),
    points0_(NullObjectRef<Field<PointType>>()),
    patch0Ptr_(nullptr),
    localPoints0Ptr_(nullptr)
{}


template<class FaceList, class PointField>
Foam::PrimitiveOldTimePatch<FaceList, PointField>::PrimitiveOldTimePatch
(
    const PrimitiveOldTimePatch<FaceList, PointField>& patch
)
:
    PrimitivePatch<FaceList, PointField>(patch),
    points0_(patch.points0_),
    patch0Ptr_(patch.patch0Ptr_, false),
    localPoints0Ptr_(nullptr)
{}


template<class FaceList, class PointField>
Foam::PrimitiveOldTimePatch<FaceList, PointField>::PrimitiveOldTimePatch
(
    PrimitiveOldTimePatch<FaceList, PointField>&& patch
)
:
    PrimitivePatch<FaceList, PointField>(move(patch)),
    points0_(move(patch.points0_)),
    patch0Ptr_(patch.patch0Ptr_),
    localPoints0Ptr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::PrimitiveOldTimePatch<FaceList, PointField>::~PrimitiveOldTimePatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class FaceList, class PointField>
const Foam::Field
<
    typename Foam::PrimitiveOldTimePatch<FaceList, PointField>::PointType
>& Foam::PrimitiveOldTimePatch<FaceList, PointField>::localPoints0() const
{
    // !!! Cannot just call patch0Ptr_->localPoints() as this would generate
    // topology in patch0Ptr_() that is already available in the base class.
    // For now, we just duplicate the implementation in PrimitivePatch.

    if (!localPoints0Ptr_)
    {
        const labelList& meshPts = this->meshPoints();

        localPoints0Ptr_ = new Field<PointType>(meshPts.size());

        Field<PointType>& locPts = *localPoints0Ptr_;

        forAll(meshPts, pointi)
        {
            locPts[pointi] = points0_[meshPts[pointi]];
        }
    }

    // But, it would be preferable to add a method to PrimitivePatch which
    // calculates the local points given a list of points different to those
    // that are stored. Then the implementations could be shared and we could
    // do this:
    /*
    if (!localPoints0Ptr_)
    {
        localPoints0Ptr_ = this->calcLocalPoints(points0_);
    }
    */

    return *localPoints0Ptr_;
}


template<class FaceList, class PointField>
const Foam::Field
<
    typename Foam::PrimitiveOldTimePatch<FaceList, PointField>::PointType
>& Foam::PrimitiveOldTimePatch<FaceList, PointField>::faceCentres0() const
{
    return patch0Ptr_->faceCentres();
}


template<class FaceList, class PointField>
const Foam::Field
<
    typename Foam::PrimitiveOldTimePatch<FaceList, PointField>::PointType
>& Foam::PrimitiveOldTimePatch<FaceList, PointField>::faceAreas0() const
{
    return patch0Ptr_->faceAreas();
}


template<class FaceList, class PointField>
const Foam::Field
<
    typename Foam::PrimitiveOldTimePatch<FaceList, PointField>::PointType
>& Foam::PrimitiveOldTimePatch<FaceList, PointField>::faceNormals0() const
{
    return patch0Ptr_->faceNormals();
}


template<class FaceList, class PointField>
const Foam::Field
<
    typename Foam::PrimitiveOldTimePatch<FaceList, PointField>::PointType
>& Foam::PrimitiveOldTimePatch<FaceList, PointField>::pointNormals0() const
{
    // !!! See comments in localPoints0. This isn't needed for now.

    NotImplemented;
    return NullObjectRef<Field<PointType>>();
}


template<class FaceList, class PointField>
void Foam::PrimitiveOldTimePatch<FaceList, PointField>::clearOut()
{
    PrimitivePatch<FaceList, PointField>::clearOut();
    if (has0()) patch0Ptr_->clearOut();

    deleteDemandDrivenData(localPoints0Ptr_);
}


template<class FaceList, class PointField>
void Foam::PrimitiveOldTimePatch<FaceList, PointField>::clearGeom()
{
    PrimitivePatch<FaceList, PointField>::clearGeom();
    if (has0()) patch0Ptr_->clearGeom();

    deleteDemandDrivenData(localPoints0Ptr_);
}


template<class FaceList, class PointField>
void Foam::PrimitiveOldTimePatch<FaceList, PointField>::movePoints0
(
    const Field<PointType>&
)
{
    if (has0()) patch0Ptr_->clearGeom();
}


// ************************************************************************* //
