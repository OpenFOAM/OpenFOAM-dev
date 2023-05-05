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

#include "layerInfo.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::layerInfo::collide() const
{
    FatalErrorInFunction
        << "Layer extrusions collided. Check the patches/zones from which "
        << "layers are being extruded and ensure that they do not point "
        << "in opposite directions."
        << exit(FatalError);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::layerInfo::layerInfo()
:
    layer_(-labelMax),
    direction_(0),
    prevFace_(-labelMax)
{}


inline Foam::layerInfo::layerInfo(const label faceLayer, const label direction)
:
    layer_(2*faceLayer - 1),
    direction_(direction),
    prevFace_(-labelMax)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::layerInfo::faceLayer() const
{
    if ((layer_ + 1) % 2 != 0)
    {
        FatalError
            << "Face layer index requested from cell layer info"
            << exit(FatalError);
    }

    return (layer_ + 1)/2;
}


inline Foam::label Foam::layerInfo::cellLayer() const
{
    if (layer_ % 2 != 0)
    {
        FatalError
            << "Cell layer index requested from face layer info"
            << exit(FatalError);
    }

    return layer_/2;
}


template<class TrackingData>
inline bool Foam::layerInfo::valid(TrackingData& td) const
{
    return layer_ != -labelMax;
}


template<class TrackingData>
inline bool Foam::layerInfo::sameGeometry
(
    const polyMesh&,
    const layerInfo& l,
    const scalar tol,
    TrackingData& td
) const
{
    return true;
}


template<class TrackingData>
inline void Foam::layerInfo::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


template<class TrackingData>
inline bool Foam::layerInfo::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label nbrFacei,
    const layerInfo& nbrInfo,
    const scalar tol,
    TrackingData& td
)
{
    const bool o = thisCelli == mesh.faceOwner()[nbrFacei];

    if (o == (nbrInfo.direction_ < 0))
    {
        if (valid(td) && prevFace_ != nbrFacei) collide();

        layer_ = valid(td) ? -labelMax : nbrInfo.layer_ + 1;
        direction_ = 0;
        prevFace_ = nbrFacei;

        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::layerInfo::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label nbrCelli,
    const layerInfo& nbrInfo,
    const scalar tol,
    TrackingData& td
)
{
    const cell& c = mesh.cells()[nbrCelli];
    const label prevFacei = nbrInfo.prevFace_;
    const label nextFacei = c.opposingFaceLabel(prevFacei, mesh.faces());

    if (nextFacei == thisFacei)
    {
        const label direction =
            mesh.faceOwner()[thisFacei] == nbrCelli ? +1 : -1;

        if (valid(td) && direction != direction_) collide();

        layer_ = valid(td) ? -labelMax : nbrInfo.layer_ + 1;
        direction_ = direction;
        prevFace_ = -labelMax;

        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::layerInfo::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const layerInfo& nbrInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label direction = -nbrInfo.direction_;

    if (valid(td) && direction != direction_) collide();

    layer_ = valid(td) ? -labelMax : nbrInfo.layer_;
    direction_ = direction;
    prevFace_ = -labelMax;

    return true;
}


template<class TrackingData>
inline bool Foam::layerInfo::equal
(
    const layerInfo& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::layerInfo::operator==(const Foam::layerInfo& rhs) const
{
    return layer_ == rhs.layer_;
}


inline bool Foam::layerInfo::operator!=(const Foam::layerInfo& rhs) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const layerInfo& l)
{
    return os << l.layer_ << token::SPACE << l.direction_;
}


inline Foam::Istream& Foam::operator>>(Istream& is, layerInfo& l)
{
    return is >> l.layer_ >> l.direction_;
}


// ************************************************************************* //
