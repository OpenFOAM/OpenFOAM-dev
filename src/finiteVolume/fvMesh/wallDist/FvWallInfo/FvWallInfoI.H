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

#include "FvWallInfo.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class WallInfo, class Derived>
template<class TrackingData>
inline void Foam::FvWallInfoBase<WallInfo, Derived>::transform
(
    const fvPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{
    WallInfo::template type<Derived>::transform
    (
        patch.patch(),
        patchFacei,
        transform,
        td
    );
}


template<class WallInfo, class Derived>
template<class TrackingData>
inline bool Foam::FvWallInfoBase<WallInfo, Derived>::updateCell
(
    const fvMesh& mesh,
    const label thisCelli,
    const labelPair& neighbourPatchAndFacei,
    const FvWallInfoBase<WallInfo, Derived>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return
        static_cast<Derived&>(*this).update
        (
            mesh.C()[thisCelli],
            static_cast<const Derived&>(neighbourWallInfo),
            tol,
            td
        );
}


template<class WallInfo, class Derived>
template<class TrackingData>
inline bool Foam::FvWallInfoBase<WallInfo, Derived>::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const label neighbourCelli,
    const FvWallInfoBase<WallInfo, Derived>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label thisPatchi = thisPatchAndFacei.first();
    const label thisFacei = thisPatchAndFacei.second();

    return
        static_cast<Derived&>(*this).update
        (
            thisPatchi == -1
          ? mesh.Cf()[thisFacei]
          : mesh.Cf().boundaryField()[thisPatchi][thisFacei],
            static_cast<const Derived&>(neighbourWallInfo),
            tol,
            td
        );
}


template<class WallInfo, class Derived>
template<class TrackingData>
inline bool Foam::FvWallInfoBase<WallInfo, Derived>::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const FvWallInfoBase<WallInfo, Derived>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label thisPatchi = thisPatchAndFacei.first();
    const label thisFacei = thisPatchAndFacei.second();

    return
        static_cast<Derived&>(*this).update
        (
            thisPatchi == -1
          ? mesh.Cf()[thisFacei]
          : mesh.Cf().boundaryField()[thisPatchi][thisFacei],
            static_cast<const Derived&>(neighbourWallInfo),
            tol,
            td
        );
}


// ************************************************************************* //
