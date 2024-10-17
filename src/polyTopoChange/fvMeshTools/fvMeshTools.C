/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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

Foam::label Foam::fvMeshTools::addPatch
(
    fvMesh& mesh,
    const polyPatch& patch
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());

    // Check if it is already there
    label patchi = polyPatches.findIndex(patch.name());
    if (patchi != -1)
    {
        return patchi;
    }

    // Append at the end ...
    label insertPatchi = polyPatches.size();

    // ... unless there are processor patches, in which case append after the
    // last global patch
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

    mesh.addPatch(insertPatchi, patch);

    return insertPatchi;
}


void Foam::fvMeshTools::addedPatches(fvMesh& mesh)
{
    mesh.addedPatches();
}


void Foam::fvMeshTools::setPatchFields
(
    fvMesh& mesh,
    const label patchi,
    const dictionary& patchFieldDict
)
{
    #define SetPatchFieldsType(Type, FieldType, Mesh) \
        setPatchFields<FieldType<Type>>(Mesh, patchi, patchFieldDict);
    FOR_ALL_FIELD_TYPES(SetPatchFieldsType, VolField, mesh);
    FOR_ALL_FIELD_TYPES(SetPatchFieldsType, SurfaceField, mesh);
    if (mesh.foundObject<pointMesh>(pointMesh::typeName))
    {
        FOR_ALL_FIELD_TYPES
        (
            SetPatchFieldsType,
            PointField,
            pointMesh::New(mesh)
        );
    }
    #undef SetPatchFieldsType
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
    // cannot use `invert(nNewPatches, oldToNew)`

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
