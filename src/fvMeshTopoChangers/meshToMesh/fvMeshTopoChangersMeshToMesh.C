/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "fvMeshTopoChangersMeshToMesh.H"
#include "polyTopoChangeMap.H"
#include "volFields.H"
#include "surfaceInterpolate.H"
#include "pointFields.H"
#include "meshToMeshAdjustTimeStepFunctionObject.H"
#include "fvMeshToFvMesh.H"
#include "cellVolumeWeightMethod.H"
#include "surfaceToVolVelocity.H"
#include "MeshToMeshMapGeometricFields.H"
#include "polyMeshMap.H"
#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(meshToMesh, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, meshToMesh, fvMesh);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::fvMeshTopoChangers::meshToMesh::interpolateUfs()
{
    // Interpolate U to Uf

    HashTable<surfaceVectorField*> Ufs
    (
        mesh().lookupClass<surfaceVectorField>()
    );

    forAllIter(HashTable<surfaceVectorField*>, Ufs, iter)
    {
        surfaceVectorField& Uf = *iter();

        const volVectorField& U = surfaceToVolVelocity(Uf);

        if (!isNull(U))
        {
            Uf.reset(fvc::interpolate(U));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::meshToMesh::meshToMesh
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshTopoChanger(mesh),
    dict_(dict),
    cuttingPatches_(dict.lookupOrDefault("cuttingPatches", wordList::null())),
    times_(dict.lookup("times")),
    timeDelta_(dict.lookup<scalar>("timeDelta")),
    timeIndex_(-1)
{
    forAll(times_, i)
    {
        timeIndices_.insert(label((times_[i] + timeDelta_/2.0)/timeDelta_));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::meshToMesh::~meshToMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::meshToMesh::update()
{
    if (timeIndex_ == -1)
    {
        const_cast<Time&>(mesh().time()).functionObjects().append
        (
            new functionObjects::meshToMeshAdjustTimeStepFunctionObject
            (
                "meshToMeshAdjustTimeStep",
                mesh().time(),
                dict_
            )
        );
    }

    bool hasChanged = false;

    // Only refine on the first call in a time-step
    if (timeIndex_ != mesh().time().timeIndex())
    {
        timeIndex_ = mesh().time().timeIndex();
    }
    else
    {
        return hasChanged;
    }

    const scalar userTime = mesh().time().userTimeValue();

    if (timeIndices_.found((userTime + timeDelta_/2)/timeDelta_))
    {
        const word meshDir = "meshToMesh_" + mesh().time().timeName(userTime);

        Info << "Mapping to mesh " << meshDir << endl;

        hasChanged = true;

        fvMesh newMesh
        (
            IOobject
            (
                meshDir,
                mesh().time().constant(),
                mesh().time(),
                IOobject::MUST_READ
            ),
            false,
            fvMesh::stitchType::none
        );

        autoPtr<Foam::fvMeshToFvMesh> mapper;

        // Create mesh-to-mesh mapper with support for cuttingPatches
        // if specified
        if (cuttingPatches_.size())
        {
            HashSet<word> cuttingPatchTable;
            forAll(cuttingPatches_, i)
            {
                cuttingPatchTable.insert(cuttingPatches_[i]);
            }

            HashTable<word> patchMap(mesh().boundary().size());

            const polyBoundaryMesh& pbm = mesh().boundaryMesh();

            forAll(pbm, i)
            {
                if
                (
                    !cuttingPatchTable.found(pbm[i].name())
                 && !isA<processorPolyPatch>(pbm[i])
                )
                {
                    patchMap.insert(pbm[i].name(), pbm[i].name());
                }
            }

            mapper = new Foam::fvMeshToFvMesh
            (
                mesh(),
                newMesh,
                cellVolumeWeightMethod::typeName,
                patchMap,
                cuttingPatches_
            );
        }
        else
        {
            mapper = new Foam::fvMeshToFvMesh
            (
                mesh(),
                newMesh,
                cellVolumeWeightMethod::typeName
            );
        }

        mesh().reset(newMesh);

        mesh().deltaCoeffs();

        // Map all the volFields in the objectRegistry
        #define mapVolFieldType(Type, nullArg)                                 \
            MeshToMeshMapVolFields<Type>(mesh(), mapper);
        FOR_ALL_FIELD_TYPES(mapVolFieldType);

        // Set all the surfaceFields in the objectRegistry to NaN
        #define NaNSurfaceFieldType(Type, nullArg)                             \
            NaNGeometricFields                                                 \
            <Type, fvsPatchField, surfaceMesh, fvPatchFieldMapper>             \
            (mesh(), mapper);
        FOR_ALL_FIELD_TYPES(NaNSurfaceFieldType);

        // Set all the pointFields in the objectRegistry to NaN
        #define NaNPointFieldType(Type, nullArg)                               \
            NaNGeometricFields                                                 \
            <Type, pointPatchField, pointMesh, pointPatchFieldMapper>          \
            (mesh(), mapper);
        FOR_ALL_FIELD_TYPES(NaNPointFieldType);

        // Interpolate U's to Uf's
        interpolateUfs();

        polyMeshMap map(mesh(), mapper());
        mesh().mapMesh(map);
    }

    return hasChanged;
}


void Foam::fvMeshTopoChangers::meshToMesh::topoChange
(
    const polyTopoChangeMap& map
)
{}


void Foam::fvMeshTopoChangers::meshToMesh::mapMesh(const polyMeshMap& map)
{}


void Foam::fvMeshTopoChangers::meshToMesh::distribute
(
    const polyDistributionMap& map
)
{}


// ************************************************************************* //
