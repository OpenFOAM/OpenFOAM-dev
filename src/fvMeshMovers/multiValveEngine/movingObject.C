/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "multiValveEngine.H"
#include "pointConstraints.H"
#include "polyCellSet.H"
#include "syncTools.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::fvMeshMovers::multiValveEngine::movingObject::calcScale
(
    const pointMesh& pMesh,
    const scalarField& pDistMoving,
    const scalarField& pDistStatic,
    scalar dMoving,
    scalar dMax,
    scalar dStatic
)
{
    scalarField& scale = scale_.primitiveFieldRef();

    forAll(scale, pi)
    {
        // If original dMoving and dStatic regions overlay set to zero.
        // This may happen when e.g. valve-piston distance becomes small.
        if
        (
            ((pDistMoving[pi] - dMoving) <= 0) &&
            ((pDistStatic[pi] - dStatic) <= 0)
        )
        {
            dMoving = 0;
            dStatic = 0;
        }

        // - Near static patches
        if (pDistStatic[pi] - dStatic <= 0)
        {
            scale[pi] = 0;
        }
        // - Near object
        else if (pDistMoving[pi] - dMoving <= 0)
        {
            scale[pi] = 1;
        }
        // - Outside outer boundary
        else if (dMax - pDistMoving[pi] <= 0)
        {
            scale[pi] = 0;
        }
        // Linear scaling from the moving patch:
        // patch->dMoving: 1
        // dMoving->min(dMax, min(dStatic - dh, dFrozenZone)): linear 1->0
        else
        {
            const scalar d1 = pDistMoving[pi] - dMoving;
            const scalar d2 = min
            (
                dMax - pDistMoving[pi],
                pDistStatic[pi] - dStatic
            );
            scale[pi] = 1 - d1/(d1 + d2 + rootVSmall);
        }
    }

    // Convert the scale function to a cosine
    if (cosine_)
    {
        scale =
            min
            (
                max
                (
                    0.5
                  - 0.5
                   *cos(scale_.primitiveField()*constant::mathematical::pi),
                    scalar(0)
                ),
                scalar(1)
            );
    }

    pointConstraints::New(pMesh).constrain(scale_);

    if (debug)
    {
        scale_.write();
    }
}


void Foam::fvMeshMovers::multiValveEngine::movingObject::transformPoints
(
    pointField& newPoints,
    const vector& translationVector
)
{
    // Explicit filtering of displacement to avoid pointMesh looping
    if (mag(translationVector) < small)
    {
        executionCount_++;
        return;
    }

    forAll(newPoints, pi)
    {
        const point displacementVector = scale_[pi]*translationVector;
        if (mag(displacementVector) > small)
        {
            newPoints[pi] += displacementVector;
        }
    }

    executionCount_++;
}


void Foam::fvMeshMovers::multiValveEngine::movingObject::createStaticPatchSet()
{
    staticPatchSet_.clear();

    forAll(meshMover_.mesh().boundaryMesh(), patchI)
    {
        const polyPatch& pp = meshMover_.mesh().boundaryMesh()[patchI];

        // Exclude non-static patches
        if
        (
           !polyPatch::constraintType(pp.type())
        && !meshMover_.slidingPatchSet_.found(pp.index())
        && !patchSet.found(pp.index())
        )
        {
            staticPatchSet_.insert(pp.index());
        }
    }
}


void Foam::fvMeshMovers::multiValveEngine::movingObject::initPatchSets()
{
    // Set patch-sets
    patchSet_ = meshMover_.mesh().boundaryMesh().patchSet(patchNames_);
    createStaticPatchSet();
}


Foam::labelHashSet
Foam::fvMeshMovers::multiValveEngine::movingObject::movingPointZones() const
{
    labelHashSet movingPointZones;

    if (movingPointZones_.size())
    {
        forAll(movingPointZones_, i)
        {
            const labelList indices
            (
                meshMover_.mesh().pointZones().findIndices(movingPointZones_[i])
            );

            if (indices.size())
            {
                movingPointZones.insert(indices);
                Info<< "    pointZone " << movingPointZones_[i]
                    << " will move with the object " << name << endl;
            }
            else
            {
                Info<< "    movingZone " << movingPointZones_[i]
                    << " not found in pointZones" << endl;
            }
        }
    }

    return movingPointZones;
}


Foam::labelHashSet
Foam::fvMeshMovers::multiValveEngine::movingObject::staticPointZones() const
{
    labelHashSet staticPointZones;

    if (frozenPointZones_.size())
    {
        forAll(frozenPointZones_, i)
        {
            const labelList indices
            (
                meshMover_.mesh().pointZones().findIndices(frozenPointZones_[i])
            );

            if (indices.size())
            {
                staticPointZones.insert(indices);
                Info<< "    pointZone " << frozenPointZones_[i]
                    << " is frozen (stationary)" << endl;
            }
            else
            {
                Info<< "    frozenZone " << frozenPointZones_[i]
                    << " not found in pointZones" << endl;
            }
        }
    }

    if (meshMover_.frozenPointZones_.size())
    {
        forAll(meshMover_.frozenPointZones_, i)
        {
            const labelList indices
            (
                meshMover_.mesh().pointZones().findIndices
                (
                    meshMover_.frozenPointZones_[i]
                )
            );

            if (indices.size())
            {
                staticPointZones.insert(indices);
                Info<< "    pointZone " << meshMover_.frozenPointZones_[i]
                    << " is frozen (stationary)" << endl;
            }
            else
            {
                Info<< "    frozenZone " << meshMover_.frozenPointZones_[i]
                    << " not found in pointZones" << endl;
            }
        }
    }

    return staticPointZones;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::multiValveEngine::movingObject::movingObject
(
    const word& objectName,
    const multiValveEngine& engine,
    const dictionary& dict
)
:
    meshMover_(engine),
    name(objectName),
    axis(dict.lookup<vector>("axis", dimless)),
    motion_(Function1<scalar>::New("motion", unitNone, dimLength, dict)),
    patchNames_(dict.lookup("patches")),
    maxMotionDistance_
    (
        dict.lookupOrDefault<scalar>("maxMotionDistance", dimLength, great)
    ),
    movingFrozenLayerThickness_
    (
        dict.lookupOrDefault<scalar>("movingFrozenLayerThickness", dimLength, 0)
    ),
    staticFrozenLayerThickness_
    (
        dict.lookupOrDefault<scalar>("staticFrozenLayerThickness", dimLength, 0)
    ),
    movingPointZones_
    (
        dict.lookupOrDefault("movingZones", wordReList::null())
    ),
    frozenPointZones_
    (
        dict.lookupOrDefault("frozenZones", wordReList::null())
    ),
    scale_
    (
        IOobject
        (
            "motionScale_" + name,
            meshMover_.mesh().time().name(),
            meshMover_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointMesh::New(meshMover_.mesh()),
        dimensionedScalar(dimless, 0)
    ),
    cosine_(dict.lookupOrDefault("cosineScaling", false)),
    travelInterval_
    (
        dict.lookupOrDefault<scalar>("travelInterval", dimLength, great)
    ),
    executionCount_(0),
    position0_(-great),
    patchSet(patchSet_)
{
    Info << indent << "Setting motion for " << name << endl;

    initPatchSets();

    if (patchSet_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Empty patchSet" << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshMovers::multiValveEngine::movingObject::mapMesh
(
    const polyMeshMap&
)
{
    // Reset patch sets.
    initPatchSets();

    // scale_ is resized by the meshToMesh mapper
    scale_ = Zero;

    // Reset count of the scale_ field updates
    executionCount_ = 0;

    // Reset position at last scale update
    position0_ = -great;
}


// ************************************************************************* //
