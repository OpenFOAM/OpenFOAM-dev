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

#include "WallInfo.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class WallLocation>
template<class TrackingData>
inline bool Foam::WallInfo<WallLocation>::sameGeometry
(
    const polyMesh& mesh,
    const WallInfo<WallLocation>& w2,
    const scalar tol,
    TrackingData& td
) const
{
    return WallLocation::sameGeometry(w2, tol, td);
}


template<class WallLocation>
template<class TrackingData>
inline void Foam::WallInfo<WallLocation>::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{
    WallLocation::transform(transform, td);
}


template<class WallLocation>
template<class TrackingData>
inline bool Foam::WallInfo<WallLocation>::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label neighbourFacei,
    const WallInfo<WallLocation>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return
        WallLocation::update
        (
            mesh.cellCentres()[thisCelli],
            neighbourWallInfo,
            tol,
            td
        );
}


template<class WallLocation>
template<class TrackingData>
inline bool Foam::WallInfo<WallLocation>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const WallInfo<WallLocation>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return
        WallLocation::update
        (
            mesh.faceCentres()[thisFacei],
            neighbourWallInfo,
            tol,
            td
        );
}


template<class WallLocation>
template<class TrackingData>
inline bool Foam::WallInfo<WallLocation>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const WallInfo<WallLocation>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return
        WallLocation::update
        (
            mesh.faceCentres()[thisFacei],
            neighbourWallInfo,
            tol,
            td
        );
}


// ************************************************************************* //
