/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "DeltaInfoData.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::DeltaInfoData<Type>::DeltaInfoData()
:
    direction_(0),
    prevFace_(-labelMax),
    data_()
{}


template<class Type>
inline Foam::DeltaInfoData<Type>::DeltaInfoData
(
    const label direction,
    const Type& data
)
:
    direction_(direction),
    prevFace_(-labelMax),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline const Type& Foam::DeltaInfoData<Type>::data() const
{
    return data_;
}


template<class Type>
template<class TrackingData>
inline bool Foam::DeltaInfoData<Type>::valid(TrackingData& td) const
{
    return direction_ != 0 || prevFace_ != -labelMax;
}


template<class Type>
template<class TrackingData>
inline bool Foam::DeltaInfoData<Type>::sameGeometry
(
    const polyMesh&,
    const DeltaInfoData<Type>& l,
    const scalar tol,
    TrackingData& td
) const
{
    return true;
}


template<class Type>
template<class TrackingData>
inline void Foam::DeltaInfoData<Type>::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


template<class Type>
template<class TrackingData>
inline bool Foam::DeltaInfoData<Type>::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label neighbourFacei,
    const DeltaInfoData<Type>& neighbourDeltaInfo,
    const scalar tol,
    TrackingData& td
)
{
    const bool o = thisCelli == mesh.faceOwner()[neighbourFacei];

    if (o == (neighbourDeltaInfo.direction_ < 0))
    {
        direction_ = 0;
        prevFace_ = neighbourFacei;
        data_ = neighbourDeltaInfo.data_;

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
template<class TrackingData>
inline bool Foam::DeltaInfoData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const DeltaInfoData<Type>& neighbourDeltaInfo,
    const scalar tol,
    TrackingData& td
)
{
    const cell& c = mesh.cells()[neighbourCelli];
    const label prevFacei = neighbourDeltaInfo.prevFace_;
    const label nextFacei = c.opposingFaceLabel(prevFacei, mesh.faces());

    if (nextFacei == thisFacei)
    {
        const label direction =
            mesh.faceOwner()[thisFacei] == neighbourCelli ? +1 : -1;

        direction_ = direction;
        prevFace_ = -labelMax;
        data_ = neighbourDeltaInfo.data_;

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
template<class TrackingData>
inline bool Foam::DeltaInfoData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const DeltaInfoData<Type>& neighbourDeltaInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label direction = -neighbourDeltaInfo.direction_;

    direction_ = direction;
    prevFace_ = -labelMax;
    data_ = neighbourDeltaInfo.data_;

    return true;
}


template<class Type>
template<class TrackingData>
inline bool Foam::DeltaInfoData<Type>::equal
(
    const DeltaInfoData<Type>& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::DeltaInfoData<Type>::operator==
(
    const Foam::DeltaInfoData<Type>& rhs
) const
{
    return direction_ == rhs.direction_ && prevFace_ == rhs.prevFace_;
}


template<class Type>
inline bool Foam::DeltaInfoData<Type>::operator!=
(
    const Foam::DeltaInfoData<Type>& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const DeltaInfoData<Type>& l)
{
    return os
        << l.direction_ << token::SPACE
        << l.prevFace_ << token::SPACE
        << l.data_;
}


template<class Type>
Foam::Istream& Foam::operator>>(Istream& is, DeltaInfoData<Type>& l)
{
    return is >> l.direction_ >> l.prevFace_ >> l.data_;
}


// ************************************************************************* //
