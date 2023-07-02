/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
#include "intersectionCellsToCells.H"
#include "surfaceToVolVelocity.H"
#include "MeshToMeshMapGeometricFields.H"
#include "polyMeshMap.H"
#include "processorPolyPatch.H"
#include "setSizeFvPatchFieldMapper.H"
#include "setSizePointPatchFieldMapper.H"
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

bool Foam::fvMeshTopoChangers::meshToMesh::forward() const
{
    return
        cycle_ == 0
     || int((mesh().time().userTimeValue() - begin_)/cycle_) % 2 == 0;
}


Foam::scalar Foam::fvMeshTopoChangers::meshToMesh::meshTime() const
{
    const Time& time = mesh().time();

    if (repeat_ > 0)
    {
        return begin_ + fmod(time.userTimeValue() - begin_, repeat_);
    }
    else if (cycle_ > 0)
    {
        if (forward())
        {
            return
                begin_
              + fmod(time.userTimeValue() - begin_, cycle_);
        }
        else
        {
            return
                begin_ + cycle_
              - fmod(time.userTimeValue() - begin_, cycle_);
        }
    }
    else
    {
        return time.userTimeValue();
    }
}


void Foam::fvMeshTopoChangers::meshToMesh::interpolateUfs()
{
    // Interpolate U to Uf
    UPtrList<surfaceVectorField> Ufs(mesh().curFields<surfaceVectorField>());

    forAll(Ufs, i)
    {
        surfaceVectorField& Uf = Ufs[i];

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
    times_(dict.lookup("times")),
    timeDelta_(dict.lookup<scalar>("timeDelta")),
    begin_(dict.lookupOrDefault("begin", mesh().time().beginTime().value())),
    repeat_(dict.lookupOrDefault("repeat", 0.0)),
    cycle_(dict.lookupOrDefault("cycle", 0.0)),
    timeIndex_(-1)
{
    if (repeat_ > 0 && cycle_ > 0)
    {
        FatalIOErrorInFunction(dict)
            << "Both 'repeat' and 'cycle' options specified"
            << exit(FatalIOError);
    }

    forAll(times_, i)
    {
        timeIndices_.insert(int64_t((times_[i] + timeDelta_/2.0)/timeDelta_));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::meshToMesh::~meshToMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fvMeshTopoChangers::meshToMesh::timeToNextMesh() const
{
    const Time& time = mesh().time();

    if
    (
        repeat_ > 0
     || cycle_ > 0
     || time.userTimeValue() + timeDelta_ < times_.last()
    )
    {
        const scalar meshTime = this->meshTime();

        if (cycle_ == 0 || int((time.userTimeValue() - begin_)/cycle_) % 2 == 0)
        {
            forAll(times_, i)
            {
                if (times_[i] > meshTime + timeDelta_)
                {
                    return time.userTimeToTime(times_[i] - meshTime);
                }
            }
        }
        else
        {
            forAllReverse(times_, i)
            {
                if (times_[i] < meshTime - timeDelta_)
                {
                    return time.userTimeToTime(meshTime - times_[i]);
                }
            }
        }
    }

    return vGreat;
}


bool Foam::fvMeshTopoChangers::meshToMesh::update()
{
    const Time& time = mesh().time();

    // Add the meshToMeshAdjustTimeStepFunctionObject functionObject
    // if not already present
    if
    (
        time.functionObjects().findObjectID("meshToMeshAdjustTimeStep")
     == -1
    )
    {
        const_cast<Time&>(time).functionObjects().append
        (
            new functionObjects::meshToMeshAdjustTimeStepFunctionObject
            (
                "meshToMeshAdjustTimeStep",
                time,
                dict_
            )
        );
    }

    // Only refine on the first call in a time-step
    if (timeIndex_ != time.timeIndex())
    {
        timeIndex_ = time.timeIndex();
    }
    else
    {
        return false;
    }

    // Obtain the mesh time from the user time
    // repeating or cycling the mesh sequence if required

    const scalar meshTime = this->meshTime();

    if (timeIndices_.found((meshTime + timeDelta_/2)/timeDelta_))
    {
        const word otherMeshDir =
            "meshToMesh_" + time.timeName(meshTime);

        Info << "Mapping to mesh " << otherMeshDir << endl;

        fvMesh otherMesh
        (
            IOobject
            (
                otherMeshDir,
                time.constant(),
                time,
                IOobject::MUST_READ
            ),
            false,
            fvMesh::stitchType::none
        );

        mesh().swap(otherMesh);

        fvMeshToFvMesh mapper
        (
            otherMesh,
            mesh(),
            cellsToCellss::intersection::typeName
        );

        // Ensure the deltaCoeffs are available for constraint patch evaluation
        mesh().deltaCoeffs();

        // Map all the volFields in the objectRegistry
        #define mapVolFieldType(Type, nullArg)                                 \
            MeshToMeshMapVolFields<Type>(mesh(), mapper);
        FOR_ALL_FIELD_TYPES(mapVolFieldType);

        // Map all the volFields in the objectRegistry
        #define mapVolInternalFieldType(Type, nullArg)                         \
            MeshToMeshMapVolInternalFields<Type>(mesh(), mapper);
        FOR_ALL_FIELD_TYPES(mapVolInternalFieldType);

        // Set all the surfaceFields in the objectRegistry to NaN
        #define NaNSurfaceFieldType(Type, nullArg)                             \
            NaNGeometricFields                                                 \
            <Type, fvsPatchField, surfaceMesh, setSizeFvPatchFieldMapper>      \
            (mesh(), mapper);
        FOR_ALL_FIELD_TYPES(NaNSurfaceFieldType);

        // Set all the pointFields in the objectRegistry to NaN
        #define NaNPointFieldType(Type, nullArg)                               \
            NaNGeometricFields                                                 \
            <Type, pointPatchField, pointMesh, setSizePointPatchFieldMapper>   \
            (mesh(), mapper);
        FOR_ALL_FIELD_TYPES(NaNPointFieldType);

        // Interpolate U's to Uf's
        interpolateUfs();

        polyMeshMap map(mesh(), mapper);

        mesh().mapMesh(map);

        return true;
    }
    else
    {
        return false;
    }
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
