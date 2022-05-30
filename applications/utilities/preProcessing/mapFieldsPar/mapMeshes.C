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

#include "mapMeshes.H"
#include "surfaceMesh.H"
#include "processorFvPatch.H"
#include "pointMesh.H"
#include "mapLagrangian.H"
#include "MapVolFields.H"
#include "UnMapped.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::mapMesh
(
    const meshToMesh& interp,
    const HashSet<word>& selectedFields,
    const bool noLagrangian
)
{
    {
        const polyMesh& meshSource = interp.srcRegion();

        // Search for list of objects for this time
        IOobjectList objects(meshSource, meshSource.time().timeName());

        // Map volFields
        // ~~~~~~~~~~~~~
        MapVolFields<scalar>
        (
            objects,
            selectedFields,
            interp
        );

        MapVolFields<vector>
        (
            objects,
            selectedFields,
            interp
        );
        MapVolFields<sphericalTensor>
        (
            objects,
            selectedFields,
            interp
        );
        MapVolFields<symmTensor>
        (
            objects,
            selectedFields,
            interp
        );
        MapVolFields<tensor>
        (
            objects,
            selectedFields,
            interp
        );
    }

    {
        const polyMesh& meshTarget = interp.tgtRegion();

        // Search for list of target objects for this time
        IOobjectList objects(meshTarget, meshTarget.time().timeName());

        // Mark surfaceFields as unmapped
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        UnMapped<surfaceScalarField>(objects);
        UnMapped<surfaceVectorField>(objects);
        UnMapped<surfaceSphericalTensorField>(objects);
        UnMapped<surfaceSymmTensorField>(objects);
        UnMapped<surfaceTensorField>(objects);

        // Mark pointFields as unmapped
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        UnMapped<pointScalarField>(objects);
        UnMapped<pointVectorField>(objects);
        UnMapped<pointSphericalTensorField>(objects);
        UnMapped<pointSymmTensorField>(objects);
        UnMapped<pointTensorField>(objects);
    }

    if (!noLagrangian)
    {
        mapLagrangian(interp);
    }
}


void Foam::mapConsistentMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const word& mapMethod,
    const HashSet<word>& selectedFields,
    const bool noLagrangian
)
{
    Info<< nl << "Consistently creating and mapping fields for time "
        << meshSource.time().timeName() << nl << endl;

    meshToMesh interp(meshSource, meshTarget, mapMethod);

    mapMesh
    (
        interp,
        selectedFields,
        noLagrangian
    );
}


void Foam::mapSubMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches,
    const word& mapMethod,
    const HashSet<word>& selectedFields,
    const bool noLagrangian
)
{
    Info<< nl << "Creating and mapping fields for time "
        << meshSource.time().timeName() << nl << endl;

    meshToMesh interp
    (
        meshSource,
        meshTarget,
        mapMethod,
        patchMap,
        cuttingPatches
    );

    mapMesh
    (
        interp,
        selectedFields,
        noLagrangian
    );
}


// ************************************************************************* //
