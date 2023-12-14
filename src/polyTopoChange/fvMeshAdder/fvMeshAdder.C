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

#include "fvMesh.H"
#include "fvMeshAdder.H"
#include "faceCoupleInfo.H"
#include "fvMesh.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(fvMeshAdder, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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
        label newFacei = oldToNew[oldStart+i];

        if (newFacei >= newStart && newFacei < newStart+newSize)
        {
            newToOld[newFacei-newStart] = i;
        }
    }
    return newToOld;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapAddedPolyMesh> Foam::fvMeshAdder::add
(
    fvMesh& mesh0,
    const fvMesh& mesh1,
    const faceCoupleInfo& coupleInfo,
    const bool validBoundary
)
{
    // Store old mesh0 point maps
    labelListList oldMeshPoints0;
    const bool havePointMesh =
        mesh0.foundObject<pointMesh>(pointMesh::typeName);
    if (havePointMesh)
    {
        const polyBoundaryMesh& pbm0 = mesh0.boundaryMesh();
        oldMeshPoints0.setSize(pbm0.size());
        forAll(pbm0, patchi)
        {
            oldMeshPoints0[patchi] = pbm0[patchi].meshPoints();
        }
    }

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
    forAll(patches, patchi)
    {
        fvPatches.set(patchi, fvPatch::New(patches[patchi], fvPatches));
    }

    // Do the mapping of the stored fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MapVolFields<scalar>(mapPtr, mesh0, mesh1);
    MapVolFields<vector>(mapPtr, mesh0, mesh1);
    MapVolFields<sphericalTensor>(mapPtr, mesh0, mesh1);
    MapVolFields<symmTensor>(mapPtr, mesh0, mesh1);
    MapVolFields<tensor>(mapPtr, mesh0, mesh1);

    MapSurfaceFields<scalar>(mapPtr, mesh0, mesh1);
    MapSurfaceFields<vector>(mapPtr, mesh0, mesh1);
    MapSurfaceFields<sphericalTensor>(mapPtr, mesh0, mesh1);
    MapSurfaceFields<symmTensor>(mapPtr, mesh0, mesh1);
    MapSurfaceFields<tensor>(mapPtr, mesh0, mesh1);

    if (havePointMesh)
    {
        // Recreate point mesh
        const pointMesh& pointMesh0 = pointMesh::New(mesh0);

        MapPointFields<scalar>(mapPtr, pointMesh0, oldMeshPoints0, mesh1);
        MapPointFields<vector>(mapPtr, pointMesh0, oldMeshPoints0, mesh1);
        MapPointFields<sphericalTensor>
        (
            mapPtr,
            pointMesh0,
            oldMeshPoints0,
            mesh1
        );
        MapPointFields<symmTensor>(mapPtr, pointMesh0, oldMeshPoints0, mesh1);
        MapPointFields<tensor>(mapPtr, pointMesh0, oldMeshPoints0, mesh1);
    }

    MapDimFields<scalar>(mapPtr, mesh0, mesh1);
    MapDimFields<vector>(mapPtr, mesh0, mesh1);
    MapDimFields<sphericalTensor>(mapPtr, mesh0, mesh1);
    MapDimFields<symmTensor>(mapPtr, mesh0, mesh1);
    MapDimFields<tensor>(mapPtr, mesh0, mesh1);

    return mapPtr;
}


// ************************************************************************* //
