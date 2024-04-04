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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::multiValveEngine::valveObject::valveObject
(
    const word& name,
    const multiValveEngine& engine,
    const dictionary& dict
)
:
    movingObject(name, engine, dict),
    minLift_(dict.lookup<scalar>("minLift"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fvMeshMovers::multiValveEngine::valveObject::lift
(
    const scalar theta
) const
{
    return motion_->value(theta);
}


bool Foam::fvMeshMovers::multiValveEngine::valveObject::isOpen() const
{
    return lift(meshMover_.userTime()) >= minLift_;
}


Foam::scalar Foam::fvMeshMovers::multiValveEngine::valveObject::lift() const
{
    return max
    (
        lift(meshMover_.userTime()),
        minLift_
    );
}


Foam::scalar Foam::fvMeshMovers::multiValveEngine::valveObject::speed() const
{
    return
       -(
            lift()
          - max
            (
                lift(meshMover_.userTime() - meshMover_.userDeltaT()),
                minLift_
            )
        )/(meshMover_.mesh().time().deltaTValue() + vSmall);
}


Foam::scalar
Foam::fvMeshMovers::multiValveEngine::valveObject::displacement() const
{
    return
        lift(meshMover_.userTime() - meshMover_.userDeltaT())
      - lift(meshMover_.userTime());
}


void Foam::fvMeshMovers::multiValveEngine::valveObject::updatePoints
(
    pointField& newPoints
)
{
    // Update points only if valve is moving
    if (isOpen())
    {
        const scalar position = this->lift();

        // Update a cached scale_ field if needed
        if (mag(position - position0_) > travelInterval_)
        {
            const pointMesh& pMesh = pointMesh::New(meshMover_.mesh());
            const pointField& points(meshMover_.mesh().points());

            const pointDist pDistMoving
            (
                pMesh,
                patchSet,
                movingPointZones(),
                points,
                maxMotionDistance_
            );

            const pointDist pDistStatic
            (
                pMesh,
                staticPatchSet_,
                staticPointZones(),
                points,
                maxMotionDistance_
            );

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
}


// ************************************************************************* //
