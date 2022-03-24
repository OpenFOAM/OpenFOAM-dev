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

#include "patchDistFuncs.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<template<class> class MapType>
void Foam::patchDistFuncs::correctBoundaryFaceCells
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& wallDistCorrected,
    MapType<label>& nearestFace
)
{
    // Size neighbours array for maximum possible (= size of largest patch)
    DynamicList<label> neighbours(maxPatchSize(mesh, patchIDs));

    // Correct all cells with face on wall
    const vectorField& cellCentres = mesh.cellCentres();

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const label patchi = iter.key();
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        // Check cells with face on wall
        forAll(patch, patchFacei)
        {
            const label celli = patch.faceCells()[patchFacei];

            getPointNeighbours
            (
                patch,
                patchFacei,
                neighbours
            );

            label minFacei = -1;

            wallDistCorrected[celli] =
                smallestDist
                (
                    cellCentres[celli],
                    patch,
                    neighbours,
                    minFacei
                );

            // Store wallCell and its nearest neighbour
            nearestFace.insert(celli, minFacei);
        }
    }
}


template<template<class> class PatchField, template<class> class MapType>
void Foam::patchDistFuncs::correctBoundaryFaceFaceCells
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    FieldField<PatchField, scalar>& wallDistCorrected,
    MapType<label>& nearestFace
)
{
    // Size neighbours array for maximum possible (= size of largest patch)
    DynamicList<label> neighbours(maxPatchSize(mesh, patchIDs));

    // Correct all faces on a wall
    const vectorField& cellCentres = mesh.cellCentres();

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const label patchi = iter.key();
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        // Check cells with face on wall
        forAll(patch, patchFacei)
        {
            const label celli = patch.faceCells()[patchFacei];

            getPointNeighbours
            (
                patch,
                patchFacei,
                neighbours
            );

            label minFacei = -1;

            wallDistCorrected[patchi][patchFacei] =
                smallestDist
                (
                    cellCentres[celli],
                    patch,
                    neighbours,
                    minFacei
                );

            // Store wallCell and its nearest neighbour
            nearestFace.insert(celli, minFacei);
        }
    }
}


// ************************************************************************* //
