/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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
#include "processorPolyPatch.H"
#include "pointFields.H"

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

    if (!isA<processorPolyPatch>(patch))
    {
        forAll(polyPatches, patchi)
        {
            const polyPatch& pp = polyPatches[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                insertPatchi = patchi;
                break;
            }
        }
    }

    mesh.addPatch
    (
        insertPatchi,
        patch,
        patchFieldDict,
        defaultPatchFieldType,
        validBoundary
    );

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
    setPatchFields<surfaceSphericalTensorField>(mesh, patchi, patchFieldDict);
    setPatchFields<surfaceSymmTensorField>(mesh, patchi, patchFieldDict);
    setPatchFields<surfaceTensorField>(mesh, patchi, patchFieldDict);

    if (mesh.foundObject<pointMesh>(pointMesh::typeName))
    {
        pointMesh& pm = const_cast<pointMesh&>(pointMesh::New(mesh));
        setPatchFields<pointScalarField>(pm, patchi, patchFieldDict);
        setPatchFields<pointVectorField>(pm, patchi, patchFieldDict);
        setPatchFields<pointSphericalTensorField>(pm, patchi, patchFieldDict);
        setPatchFields<pointSymmTensorField>(pm, patchi, patchFieldDict);
        setPatchFields<pointTensorField>(pm, patchi, patchFieldDict);
    }
}


void Foam::fvMeshTools::zeroPatchFields(fvMesh& mesh, const label patchi)
{
    setPatchFields<volScalarField>(mesh, patchi, Zero);
    setPatchFields<volVectorField>(mesh, patchi, Zero);
    setPatchFields<volSphericalTensorField>(mesh, patchi, Zero);
    setPatchFields<volSymmTensorField>(mesh, patchi, Zero);
    setPatchFields<volTensorField>(mesh, patchi, Zero);

    setPatchFields<surfaceScalarField>(mesh, patchi, Zero);
    setPatchFields<surfaceVectorField>(mesh, patchi, Zero);
    setPatchFields<surfaceSphericalTensorField>(mesh, patchi, Zero);
    setPatchFields<surfaceSymmTensorField>(mesh, patchi, Zero);
    setPatchFields<surfaceTensorField>(mesh, patchi, Zero);

    if (mesh.foundObject<pointMesh>(pointMesh::typeName))
    {
        pointMesh& pm = const_cast<pointMesh&>(pointMesh::New(mesh));
        setPatchFields<pointScalarField>(pm, patchi, Zero);
        setPatchFields<pointVectorField>(pm, patchi, Zero);
        setPatchFields<pointSphericalTensorField>(pm, patchi, Zero);
        setPatchFields<pointSymmTensorField>(pm, patchi, Zero);
        setPatchFields<pointTensorField>(pm, patchi, Zero);
    }
}


void Foam::fvMeshTools::reorderPatches
(
    fvMesh& mesh,
    const labelList& oldToNew,
    const label nNewPatches,
    const bool validBoundary
)
{
    // Note: oldToNew might have entries beyond nNewPatches so
    // cannot use invert
    //const labelList newToOld(invert(nNewPatches, oldToNew));
    labelList newToOld(nNewPatches, -1);
    forAll(oldToNew, i)
    {
        label newi = oldToNew[i];

        if (newi >= 0 && newi < nNewPatches)
        {
            newToOld[newi] = i;
        }
    }
    mesh.reorderPatches(newToOld, validBoundary);
}


// ************************************************************************* //
