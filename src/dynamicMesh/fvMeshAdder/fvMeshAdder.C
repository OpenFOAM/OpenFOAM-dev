/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "fvMesh.H"
#include "fvMeshAdder.H"
#include "faceCoupleInfo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Calculate map from new patch faces to old patch faces. -1 where
//  could not map.
Foam::labelList Foam::fvMeshAdder::calcPatchMap
(
    const label oldStart,
    const label oldSize,
    const labelList& oldToNew,
    const polyPatch& newPatch,
    const label unmappedValue
)
{
    labelList newToOld(newPatch.size(), unmappedValue);

    label newStart = newPatch.start();
    label newSize = newPatch.size();

    for (label i = 0; i < oldSize; i++)
    {
        label newFaceI = oldToNew[oldStart+i];

        if (newFaceI >= newStart && newFaceI < newStart+newSize)
        {
            newToOld[newFaceI-newStart] = i;
        }
    }
    return newToOld;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Inplace add mesh1 to mesh0
Foam::autoPtr<Foam::mapAddedPolyMesh> Foam::fvMeshAdder::add
(
    fvMesh& mesh0,
    const fvMesh& mesh1,
    const faceCoupleInfo& coupleInfo,
    const bool validBoundary
)
{
    mesh0.clearOut();

    // Resulting merged mesh (polyMesh only!)
    autoPtr<mapAddedPolyMesh> mapPtr
    (
        polyMeshAdder::add
        (
            mesh0,
            mesh1,
            coupleInfo,
            validBoundary
        )
    );

    // Adjust the fvMesh part.
    const polyBoundaryMesh& patches = mesh0.boundaryMesh();

    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh0.boundary());
    fvPatches.setSize(patches.size());
    forAll(patches, patchI)
    {
        fvPatches.set(patchI, fvPatch::New(patches[patchI], fvPatches));
    }

    // Do the mapping of the stored fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fvMeshAdder::MapVolFields<scalar>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapVolFields<vector>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapVolFields<sphericalTensor>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapVolFields<symmTensor>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapVolFields<tensor>(mapPtr, mesh0, mesh1);

    fvMeshAdder::MapSurfaceFields<scalar>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapSurfaceFields<vector>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapSurfaceFields<sphericalTensor>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapSurfaceFields<symmTensor>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapSurfaceFields<tensor>(mapPtr, mesh0, mesh1);

    return mapPtr;
}


// ************************************************************************* //
