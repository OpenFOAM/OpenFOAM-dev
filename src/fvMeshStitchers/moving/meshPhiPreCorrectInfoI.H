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

#include "meshPhiPreCorrectInfo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::meshPhiPreCorrectInfo::meshPhiPreCorrectInfo()
:
    layer_(-labelMax),
    weight_(NaN)
{}


inline Foam::meshPhiPreCorrectInfo::meshPhiPreCorrectInfo
(
    const label layer,
    const scalar weight
)
:
    layer_(layer),
    weight_(weight)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::meshPhiPreCorrectInfo::layer() const
{
    return layer_;
}


inline Foam::scalar Foam::meshPhiPreCorrectInfo::weight() const
{
    return weight_;
}


template<class TrackingData>
inline bool Foam::meshPhiPreCorrectInfo::valid(TrackingData& td) const
{
    return layer_ != -labelMax;
}


template<class TrackingData>
inline bool Foam::meshPhiPreCorrectInfo::sameGeometry
(
    const fvMesh&,
    const meshPhiPreCorrectInfo& l,
    const scalar tol,
    TrackingData& td
) const
{
    return layer_ == l.layer_;
}


template<class TrackingData>
inline void Foam::meshPhiPreCorrectInfo::transform
(
    const fvPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


template<class TrackingData>
inline bool Foam::meshPhiPreCorrectInfo::updateCell
(
    const fvMesh& mesh,
    const label thisCelli,
    const labelPair& neighbourPatchAndFacei,
    const meshPhiPreCorrectInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        layer_ = neighbourInfo.layer_ + 1;
        weight_ = 0;
    }

    if (layer_ > neighbourInfo.layer_)
    {
        weight_ += neighbourInfo.weight_;

        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::meshPhiPreCorrectInfo::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const label neighbourCelli,
    const meshPhiPreCorrectInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        layer_ = neighbourInfo.layer_ + 1;
    }

    if (layer_ > neighbourInfo.layer_)
    {
        weight_ = neighbourInfo.weight_;

        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::meshPhiPreCorrectInfo::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const meshPhiPreCorrectInfo& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (layer_ < neighbourInfo.layer_)
    {
        layer_ = neighbourInfo.layer_;
        weight_ = neighbourInfo.weight_;
    }

    return true;
}


template<class TrackingData>
inline bool Foam::meshPhiPreCorrectInfo::equal
(
    const meshPhiPreCorrectInfo& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::meshPhiPreCorrectInfo::operator==
(
    const Foam::meshPhiPreCorrectInfo& rhs
) const
{
    return layer_ == rhs.layer_ && weight_ == rhs.weight_;
}


inline bool Foam::meshPhiPreCorrectInfo::operator!=
(
    const Foam::meshPhiPreCorrectInfo& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const meshPhiPreCorrectInfo& l)
{
    return os << l.layer_ << token::SPACE << l.weight_;
}


Foam::Istream& Foam::operator>>(Istream& is, meshPhiPreCorrectInfo& l)
{
    return is >> l.layer_ >> l.weight_;
}


// ************************************************************************* //
