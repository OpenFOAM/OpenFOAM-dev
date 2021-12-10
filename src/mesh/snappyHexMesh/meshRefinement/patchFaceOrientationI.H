/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "orientedSurface.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::patchFaceOrientation::patchFaceOrientation()
:
    flipStatus_(orientedSurface::UNVISITED)
{}


inline Foam::patchFaceOrientation::patchFaceOrientation
(
    const label flipStatus
)
:
    flipStatus_(flipStatus)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::patchFaceOrientation::flipStatus() const
{
    return flipStatus_;
}


inline void Foam::patchFaceOrientation::flip()
{
    if (flipStatus_ == orientedSurface::NOFLIP)
    {
        flipStatus_ = orientedSurface::FLIP;
    }
    else if (flipStatus_ == orientedSurface::FLIP)
    {
        flipStatus_ = orientedSurface::NOFLIP;
    }
}


template<class TrackingData>
inline bool Foam::patchFaceOrientation::valid(TrackingData& td) const
{
    return flipStatus_ != orientedSurface::UNVISITED;
}


template<class TrackingData>
inline void Foam::patchFaceOrientation::transform
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const tensor& rotTensor,
    const scalar tol,
    TrackingData& td
)
{}


template<class TrackingData>
inline bool Foam::patchFaceOrientation::updateEdge
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const label edgei,
    const label facei,
    const patchFaceOrientation& faceInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (valid(td))
    {
        return false;
    }

    const face& f = patch.localFaces()[facei];
    const edge& e = patch.edges()[edgei];

    patchFaceOrientation consistentInfo(faceInfo);

    // Check how edge relates to face
    if (f.edgeDirection(e) < 0)
    {
        // Create flipped version of faceInfo
        consistentInfo.flip();
    }

    operator=(consistentInfo);

    return true;
}


template<class TrackingData>
inline bool Foam::patchFaceOrientation::updateEdge
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const patchFaceOrientation& edgeinfo,
    const bool sameOrientation,
    const scalar tol,
    TrackingData& td
)
{
    if (valid(td))
    {
        return false;
    }

    // Create (flipped/unflipped) version of edgeinfo
    patchFaceOrientation consistentInfo(edgeinfo);

    if (!sameOrientation)
    {
        consistentInfo.flip();
    }

    operator=(consistentInfo);

    return true;
}


template<class TrackingData>
inline bool Foam::patchFaceOrientation::updateFace
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const label facei,
    const label edgei,
    const patchFaceOrientation& edgeinfo,
    const scalar tol,
    TrackingData& td
)
{
    if (valid(td))
    {
        return false;
    }

    // Transfer flip to face
    const face& f = patch.localFaces()[facei];
    const edge& e = patch.edges()[edgei];

    // Create (flipped/unflipped) version of edgeinfo
    patchFaceOrientation consistentInfo(edgeinfo);

    if (f.edgeDirection(e) > 0)
    {
        consistentInfo.flip();
    }

    operator=(consistentInfo);

    return true;
}


template<class TrackingData>
inline bool Foam::patchFaceOrientation::equal
(
    const patchFaceOrientation& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::patchFaceOrientation::operator==
(
    const Foam::patchFaceOrientation& rhs
) const
{
    return flipStatus() == rhs.flipStatus();
}


inline bool Foam::patchFaceOrientation::operator!=
(
    const Foam::patchFaceOrientation& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
