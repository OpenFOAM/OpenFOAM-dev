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

#include "pointEdgeLayerInfo.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::pointEdgeLayerInfo::pointEdgeLayerInfo()
:
    layer_(-labelMax)
{}


inline Foam::pointEdgeLayerInfo::pointEdgeLayerInfo(const label pointLayer)
:
    layer_(2*pointLayer - 1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::pointEdgeLayerInfo::pointLayer() const
{
    if ((layer_ + 1) % 2 != 0)
    {
        FatalError
            << "Point layer index requested from edge layer info"
            << exit(FatalError);
    }

    return (layer_ + 1)/2;
}


inline Foam::label Foam::pointEdgeLayerInfo::edgeLayer() const
{
    if (layer_ % 2 != 0)
    {
        FatalError
            << "Edge layer index requested from point layer info"
            << exit(FatalError);
    }

    return layer_/2;
}


template<class TrackingData>
inline bool Foam::pointEdgeLayerInfo::valid(TrackingData& td) const
{
    return layer_ != -labelMax;
}


template<class TrackingData>
inline bool Foam::pointEdgeLayerInfo::sameGeometry
(
    const pointEdgeLayerInfo& l,
    const scalar tol,
    TrackingData& td
) const
{
    return true;
}


template<class TrackingData>
inline void Foam::pointEdgeLayerInfo::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


template<class TrackingData>
inline bool Foam::pointEdgeLayerInfo::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgei,
    const pointEdgeLayerInfo& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        layer_ = edgeInfo.layer_ + 1;
        return true;
    }

    if (!edgeInfo.valid(td) || layer_ < edgeInfo.layer_)
    {
        return false;
    }

    FatalErrorInFunction
        << "Layer information collided. This is not a layered mesh"
        << exit(FatalError);

    return false;
}


template<class TrackingData>
inline bool Foam::pointEdgeLayerInfo::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const pointEdgeLayerInfo& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        layer_ = newPointInfo.layer_;
        return true;
    }

    FatalErrorInFunction
        << "Layer information collided. This is not a layered mesh"
        << exit(FatalError);

    return false;
}


template<class TrackingData>
inline bool Foam::pointEdgeLayerInfo::updatePoint
(
    const pointEdgeLayerInfo& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    layer_ = newPointInfo.layer_;
    return true;
}


template<class TrackingData>
inline bool Foam::pointEdgeLayerInfo::updateEdge
(
    const polyMesh& mesh,
    const label edgei,
    const label pointi,
    const pointEdgeLayerInfo& pointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        layer_ = pointInfo.layer_ + 1;
        return true;
    }

    layer_ = -labelMax;
    return false;
}


template<class TrackingData>
inline bool Foam::pointEdgeLayerInfo::equal
(
    const pointEdgeLayerInfo& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pointEdgeLayerInfo::operator==
(
    const Foam::pointEdgeLayerInfo& rhs
) const
{
    return layer_ == rhs.layer_;
}


inline bool Foam::pointEdgeLayerInfo::operator!=
(
    const Foam::pointEdgeLayerInfo& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const pointEdgeLayerInfo& l)
{
    return os << l.layer_;
}


inline Foam::Istream& Foam::operator>>(Istream& is, pointEdgeLayerInfo& l)
{
    return is >> l.layer_;
}


// ************************************************************************* //
