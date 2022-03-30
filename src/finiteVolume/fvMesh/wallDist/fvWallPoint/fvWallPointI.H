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

#include "fvWallPoint.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline void Foam::fvWallPoint::transform
(
    const fvPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{
    wallPoint::transform(patch.patch(), patchFacei, transform, td);
}


template<class TrackingData>
inline bool Foam::fvWallPoint::updateCell
(
    const fvMesh& mesh,
    const label thisCelli,
    const labelPair& neighbourPatchAndFacei,
    const fvWallPoint& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return
        this->update
        (
            mesh.C()[thisCelli],
            neighbourWallInfo,
            tol,
            td
        );
}


template<class TrackingData>
inline bool Foam::fvWallPoint::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const label neighbourCelli,
    const fvWallPoint& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label thisPatchi = thisPatchAndFacei.first();
    const label thisFacei = thisPatchAndFacei.second();

    return
        this->update
        (
            thisPatchi == -1
          ? mesh.Cf()[thisFacei]
          : mesh.Cf().boundaryField()[thisPatchi][thisFacei],
            neighbourWallInfo,
            tol,
            td
        );
}


template<class TrackingData>
inline bool Foam::fvWallPoint::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const fvWallPoint& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label thisPatchi = thisPatchAndFacei.first();
    const label thisFacei = thisPatchAndFacei.second();

    return
        this->update
        (
            thisPatchi == -1
          ? mesh.Cf()[thisFacei]
          : mesh.Cf().boundaryField()[thisPatchi][thisFacei],
            neighbourWallInfo,
            tol,
            td
        );
}


// ************************************************************************* //
