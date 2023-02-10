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

#include "patchDistWave.H"
#include "FaceCellWave.H"
#include "wallPoint.H"
#include "WallInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::labelList Foam::patchDistWave::getChangedFaces
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs
)
{
    label nChangedFaces = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        nChangedFaces += mesh.boundaryMesh()[iter.key()].size();
    }

    labelList changedFaces(nChangedFaces);
    label changedFacei = 0;

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const label patchi = iter.key();

        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        forAll(patch.faceCentres(), patchFacei)
        {
            const label meshFacei = patch.start() + patchFacei;

            changedFaces[changedFacei] = meshFacei;

            changedFacei++;
        }
    }

    return changedFaces;
}


Foam::label Foam::patchDistWave::wave
(
    const polyMesh& mesh,
    const labelList& changedFaces,
    scalarField& cellDistance
)
{
    // Initialise changedFacesInfo to face centres on patches
    List<WallInfo<wallPoint>> changedFacesInfo(changedFaces.size());
    forAll(changedFaces, changedFacei)
    {
        const label facei = changedFaces[changedFacei];

        changedFacesInfo[changedFacei] =
            WallInfo<wallPoint>(mesh.faceCentres()[facei], scalar(0));
    }

    // Do calculate patch distance by 'growing' from faces.
    List<WallInfo<wallPoint>> faceInfo(mesh.nFaces()), cellInfo(mesh.nCells());
    FaceCellWave<WallInfo<wallPoint>> wave
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        faceInfo,
        cellInfo,
        mesh.globalData().nTotalCells() + 1 // max iterations
    );

    // Copy distances into field
    label nUnset = 0;
    forAll(cellInfo, celli)
    {
        nUnset += !cellInfo[celli].valid(wave.data());

        cellDistance[celli] = cellInfo[celli].dist(wave.data());
    }

    return nUnset;
}


Foam::label Foam::patchDistWave::calculate
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& cellDistance
)
{
    return wave(mesh, getChangedFaces(mesh, patchIDs), cellDistance);
}


// ************************************************************************* //
