/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "fvMeshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Adds patch if not yet there. Returns patchID.
Foam::label Foam::fvMeshTools::addPatch
(
    fvMesh& mesh,
    const polyPatch& patch,
    const dictionary& patchFieldDict,
    const word& defaultPatchFieldType,
    const bool validBoundary
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());

    label patchi = polyPatches.findPatchID(patch.name());
    if (patchi != -1)
    {
        // Already there
        return patchi;
    }


    // Append at end unless there are processor patches
    label insertPatchi = polyPatches.size();
    label startFacei = mesh.nFaces();

    if (!isA<processorPolyPatch>(patch))
    {
        forAll(polyPatches, patchi)
        {
            const polyPatch& pp = polyPatches[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                insertPatchi = patchi;
                startFacei = pp.start();
                break;
            }
        }
    }


    // Below is all quite a hack. Feel free to change once there is a better
    // mechanism to insert and reorder patches.

    // Clear local fields and e.g. polyMesh parallelInfo.
    mesh.clearOut();

    label sz = polyPatches.size();

    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    // Add polyPatch at the end
    polyPatches.setSize(sz+1);
    polyPatches.set
    (
        sz,
        patch.clone
        (
            polyPatches,
            insertPatchi,   // index
            0,              // size
            startFacei      // start
        )
    );
    fvPatches.setSize(sz+1);
    fvPatches.set
    (
        sz,
        fvPatch::New
        (
            polyPatches[sz],  // point to newly added polyPatch
            mesh.boundary()
        )
    );

    addPatchFields<volScalarField>
    (
        mesh,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    addPatchFields<volVectorField>
    (
        mesh,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    addPatchFields<volSphericalTensorField>
    (
        mesh,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    addPatchFields<volSymmTensorField>
    (
        mesh,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    addPatchFields<volTensorField>
    (
        mesh,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );

    // Surface fields

    addPatchFields<surfaceScalarField>
    (
        mesh,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    addPatchFields<surfaceVectorField>
    (
        mesh,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    addPatchFields<surfaceSphericalTensorField>
    (
        mesh,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    addPatchFields<surfaceSymmTensorField>
    (
        mesh,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    addPatchFields<surfaceTensorField>
    (
        mesh,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );

    // Create reordering list
    // patches before insert position stay as is
    labelList oldToNew(sz+1);
    for (label i = 0; i < insertPatchi; i++)
    {
        oldToNew[i] = i;
    }
    // patches after insert position move one up
    for (label i = insertPatchi; i < sz; i++)
    {
        oldToNew[i] = i+1;
    }
    // appended patch gets moved to insert position
    oldToNew[sz] = insertPatchi;

    // Shuffle into place
    polyPatches.reorder(oldToNew, validBoundary);
    fvPatches.reorder(oldToNew);

    reorderPatchFields<volScalarField>(mesh, oldToNew);
    reorderPatchFields<volVectorField>(mesh, oldToNew);
    reorderPatchFields<volSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<volSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<volTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceScalarField>(mesh, oldToNew);
    reorderPatchFields<surfaceVectorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceTensorField>(mesh, oldToNew);

    return insertPatchi;
}


void Foam::fvMeshTools::setPatchFields
(
    fvMesh& mesh,
    const label patchi,
    const dictionary& patchFieldDict
)
{
    setPatchFields<volScalarField>(mesh, patchi, patchFieldDict);
    setPatchFields<volVectorField>(mesh, patchi, patchFieldDict);
    setPatchFields<volSphericalTensorField>(mesh, patchi, patchFieldDict);
    setPatchFields<volSymmTensorField>(mesh, patchi, patchFieldDict);
    setPatchFields<volTensorField>(mesh, patchi, patchFieldDict);
    setPatchFields<surfaceScalarField>(mesh, patchi, patchFieldDict);
    setPatchFields<surfaceVectorField>(mesh, patchi, patchFieldDict);
    setPatchFields<surfaceSphericalTensorField>
    (
        mesh,
        patchi,
        patchFieldDict
    );
    setPatchFields<surfaceSymmTensorField>(mesh, patchi, patchFieldDict);
    setPatchFields<surfaceTensorField>(mesh, patchi, patchFieldDict);
}


void Foam::fvMeshTools::zeroPatchFields(fvMesh& mesh, const label patchi)
{
    setPatchFields<volScalarField>(mesh, patchi, Zero);
    setPatchFields<volVectorField>(mesh, patchi, Zero);
    setPatchFields<volSphericalTensorField>
    (
        mesh,
        patchi,
        Zero
    );
    setPatchFields<volSymmTensorField>
    (
        mesh,
        patchi,
        Zero
    );
    setPatchFields<volTensorField>(mesh, patchi, Zero);
    setPatchFields<surfaceScalarField>(mesh, patchi, Zero);
    setPatchFields<surfaceVectorField>(mesh, patchi, Zero);
    setPatchFields<surfaceSphericalTensorField>
    (
        mesh,
        patchi,
        Zero
    );
    setPatchFields<surfaceSymmTensorField>
    (
        mesh,
        patchi,
        Zero
    );
    setPatchFields<surfaceTensorField>(mesh, patchi, Zero);
}


// Deletes last patch
void Foam::fvMeshTools::trimPatches(fvMesh& mesh, const label nPatches)
{
    // Clear local fields and e.g. polyMesh globalMeshData.
    mesh.clearOut();

    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    if (polyPatches.empty())
    {
        FatalErrorInFunction
            << "No patches in mesh"
            << abort(FatalError);
    }

    label nFaces = 0;
    for (label patchi = nPatches; patchi < polyPatches.size(); patchi++)
    {
        nFaces += polyPatches[patchi].size();
    }
    reduce(nFaces, sumOp<label>());

    if (nFaces)
    {
        FatalErrorInFunction
            << "There are still " << nFaces
            << " faces in " << polyPatches.size()-nPatches
            << " patches to be deleted" << abort(FatalError);
    }

    // Remove actual patches
    polyPatches.setSize(nPatches);
    fvPatches.setSize(nPatches);

    trimPatchFields<volScalarField>(mesh, nPatches);
    trimPatchFields<volVectorField>(mesh, nPatches);
    trimPatchFields<volSphericalTensorField>(mesh, nPatches);
    trimPatchFields<volSymmTensorField>(mesh, nPatches);
    trimPatchFields<volTensorField>(mesh, nPatches);

    trimPatchFields<surfaceScalarField>(mesh, nPatches);
    trimPatchFields<surfaceVectorField>(mesh, nPatches);
    trimPatchFields<surfaceSphericalTensorField>(mesh, nPatches);
    trimPatchFields<surfaceSymmTensorField>(mesh, nPatches);
    trimPatchFields<surfaceTensorField>(mesh, nPatches);
}


void Foam::fvMeshTools::reorderPatches
(
    fvMesh& mesh,
    const labelList& oldToNew,
    const label nNewPatches,
    const bool validBoundary
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    // Shuffle into place
    polyPatches.reorder(oldToNew, validBoundary);
    fvPatches.reorder(oldToNew);

    reorderPatchFields<volScalarField>(mesh, oldToNew);
    reorderPatchFields<volVectorField>(mesh, oldToNew);
    reorderPatchFields<volSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<volSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<volTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceScalarField>(mesh, oldToNew);
    reorderPatchFields<surfaceVectorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceTensorField>(mesh, oldToNew);

    // Remove last.
    trimPatches(mesh, nNewPatches);
}


// ************************************************************************* //
