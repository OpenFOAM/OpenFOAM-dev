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

#include "LayerInfoData.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::LayerInfoData<Type>::LayerInfoData()
:
    layerInfo(),
    data_()
{}


template<class Type>
inline Foam::LayerInfoData<Type>::LayerInfoData
(
    const label faceLayer,
    const label direction,
    const Type& data
)
:
    layerInfo(faceLayer, direction),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline const Type& Foam::LayerInfoData<Type>::data() const
{
    return data_;
}


template<class Type>
template<class TrackingData>
inline void Foam::LayerInfoData<Type>::transform
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
inline bool Foam::LayerInfoData<Type>::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label nbrFacei,
    const LayerInfoData<Type>& nbrInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (layerInfo::updateCell(mesh, thisCelli, nbrFacei, nbrInfo, tol, td))
    {
        data_ = nbrInfo.data_;
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
template<class TrackingData>
inline bool Foam::LayerInfoData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label nbrCelli,
    const LayerInfoData<Type>& nbrInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (layerInfo::updateFace(mesh, thisFacei, nbrCelli, nbrInfo, tol, td))
    {
        data_ = nbrInfo.data_;
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
template<class TrackingData>
inline bool Foam::LayerInfoData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const LayerInfoData<Type>& nbrInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (layerInfo::updateFace(mesh, thisFacei, nbrInfo, tol, td))
    {
        data_ = nbrInfo.data_;
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
template<class TrackingData>
inline bool Foam::LayerInfoData<Type>::equal
(
    const LayerInfoData<Type>& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::LayerInfoData<Type>::operator==
(
    const Foam::LayerInfoData<Type>& rhs
) const
{
    return layerInfo::operator==(rhs) && data_ == rhs.data_;
}


template<class Type>
inline bool Foam::LayerInfoData<Type>::operator!=
(
    const Foam::LayerInfoData<Type>& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LayerInfoData<Type>& l
)
{
    return os << static_cast<const layerInfo&>(l) << token::SPACE << l.data_;
}


template<class Type>
Foam::Istream& Foam::operator>>(Istream& is, LayerInfoData<Type>& l)
{
    return is >> static_cast<layerInfo&>(l) >> l.data_;
}


// ************************************************************************* //
