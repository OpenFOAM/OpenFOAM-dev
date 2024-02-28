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
#include "pointDist.H"

/* * * * * * * * * * * * * Static Private Member Data  * * * * * * * * * * * */

Foam::word Foam::fvMeshMovers::multiValveEngine::pistonObject::pistonBowlName
(
    "pistonBowl"
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::fvMeshMovers::multiValveEngine::pistonObject::bore() const
{
    const polyBoundaryMesh& pbm = meshMover_.mesh().boundaryMesh();

    // Find the maximum and minimum coordinates of the piston patch-set
    vector pistonMax(vector::min);
    vector pistonMin(vector::max);

    forAllConstIter(labelHashSet, patchSet, iter)
    {
        const label patchi = iter.key();
        if (pbm[patchi].localPoints().size())
        {
            pistonMax = max(pistonMax, max(pbm[patchi].localPoints()));
            pistonMin = min(pistonMin, min(pbm[patchi].localPoints()));
        }
    }

    reduce(pistonMax, maxOp<point>());
    reduce(pistonMin, minOp<point>());

    // Assuming the piston moves in the positive axis direction
    // remove the axis_ component to find the lateral extent of the piston
    return mag
    (
        (pistonMax - (axis& pistonMax)*pistonMax)
      - (pistonMin - (axis& pistonMin)*pistonMin)
    )/sqrt(2.0);
}


void Foam::fvMeshMovers::multiValveEngine::pistonObject::correctClearance
(
    pointDist& pDist
)
{
    clearance_ = great;

    forAllConstIter(labelHashSet, staticPatchSet_, iter)
    {
        const polyPatch& pp = meshMover_.mesh().boundaryMesh()[iter.key()];
        const labelList& meshPoints = pp.meshPoints();

        forAll(meshPoints, pointi)
        {
            clearance_ = min(clearance_, pDist[meshPoints[pointi]]);
        }
    }

    reduce(clearance_, minOp<scalar>());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::multiValveEngine::pistonObject::pistonObject
(
    const word& name,
    const multiValveEngine& engine,
    const dictionary& dict
)
:
    movingObject(name, engine, dict),
    bore_(bore()),
    clearance_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fvMeshMovers::multiValveEngine::pistonObject::position
(
    const scalar theta
) const
{
    return motion_->value(theta);
}


Foam::scalar
Foam::fvMeshMovers::multiValveEngine::pistonObject::position() const
{
    return position(meshMover_.userTime());
}


Foam::scalar
Foam::fvMeshMovers::multiValveEngine::pistonObject::displacement() const
{
    return
        position(meshMover_.userTime() - meshMover_.userDeltaT())
      - position();
}


Foam::scalar Foam::fvMeshMovers::multiValveEngine::pistonObject::speed() const
{
    return displacement()/(meshMover_.mesh().time().deltaTValue() + vSmall);
}


Foam::scalar
Foam::fvMeshMovers::multiValveEngine::pistonObject::clearance() const
{
    if (mag(position() - position0_)/travel_ > fractionalTravelInterval_)
    {
        return clearance_;
    }
    else
    {
        // Note, valve movement is not considered as valveSet may include
        // other valves than ones related to ports. Furthermore, this value
        // is only an estimate and updated rather frequently anyway.
        return clearance_ - displacement();
    }
}


void Foam::fvMeshMovers::multiValveEngine::pistonObject::updatePoints
(
    pointField& newPoints
)
{
    const scalar position = this->position();

    // Update a cached scale_ field if needed
    if (mag(position - position0_)/travel_ > fractionalTravelInterval_)
    {
        Info << "    Updating scale field" << endl;

        const pointMesh& pMesh = pointMesh::New(meshMover_.mesh());
        const pointField& points(meshMover_.mesh().points());

        pointDist pDistMoving
        (
            pMesh,
            patchSet,
            movingPointZones(),
            points,
            maxMotionDistance_
        );

        pointDist pDistStatic
        (
            pMesh,
            staticPatchSet_,
            staticPointZones(),
            points,
            maxMotionDistance_
        );

        // Update the clearance from the distance to piston field
        correctClearance(pDistMoving);

        calcScale
        (
            pMesh,
            pDistMoving,
            pDistStatic,
            movingFrozenLayerThickness_,
            maxMotionDistance_,
            staticFrozenLayerThickness_
        );

        position0_ = position;
    }

    const vector translationVector(displacement()*axis);
    transformPoints(newPoints, translationVector);
}


void Foam::fvMeshMovers::multiValveEngine::pistonObject::mapMesh
(
    const polyMeshMap& map
)
{
    movingObject::mapMesh(map);
}


// ************************************************************************* //
