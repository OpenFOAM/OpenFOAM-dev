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

#include "engineValve.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::engineValve::adjustCrankAngle(const scalar theta) const
{
    if (theta < liftProfileStart_)
    {
        scalar adjustedTheta = theta;

        while (adjustedTheta < liftProfileStart_)
        {
            adjustedTheta += liftProfileEnd_ - liftProfileStart_;
        }

        return adjustedTheta;
    }
    else if (theta > liftProfileEnd_)
    {
        scalar adjustedTheta = theta;

        while (adjustedTheta > liftProfileEnd_)
        {
            adjustedTheta -= liftProfileEnd_ - liftProfileStart_;
        }

        return adjustedTheta;
    }
    else
    {
        return theta;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::engineValve::engineValve
(
    const word& name,
    const fvMeshMover& meshMover,
    const autoPtr<coordinateSystem>& valveCS,
    const word& bottomPatchName,
    const word& poppetPatchName,
    const word& stemPatchName,
    const word& curtainInPortPatchName,
    const word& curtainInCylinderPatchName,
    const word& detachInCylinderPatchName,
    const word& detachInPortPatchName,
    const labelList& detachFaces,
    const Function1s::Table<scalar>& liftProfile,
    const scalar minLift,
    const scalar minTopLayer,
    const scalar maxTopLayer,
    const scalar minBottomLayer,
    const scalar maxBottomLayer,
    const scalar diameter
)
:
    name_(name),
    meshMover_(refCast<const fvMeshMovers::engine>(meshMover)),
    csPtr_(valveCS),
    bottomPatch_(bottomPatchName, meshMover_.mesh().boundaryMesh()),
    poppetPatch_(poppetPatchName, meshMover_.mesh().boundaryMesh()),
    stemPatch_(stemPatchName, meshMover_.mesh().boundaryMesh()),
    curtainInPortPatch_
    (
        curtainInPortPatchName, meshMover_.mesh().boundaryMesh()
    ),
    curtainInCylinderPatch_
    (
        curtainInCylinderPatchName, meshMover_.mesh().boundaryMesh()
    ),
    detachInCylinderPatch_
    (
        detachInCylinderPatchName, meshMover_.mesh().boundaryMesh()
    ),
    detachInPortPatch_(detachInPortPatchName, meshMover_.mesh().boundaryMesh()),
    detachFaces_(detachFaces),
    liftProfile_(liftProfile),
    liftProfileStart_(min(liftProfile_.x())),
    liftProfileEnd_(max(liftProfile_.x())),
    minLift_(minLift),
    minTopLayer_(minTopLayer),
    maxTopLayer_(maxTopLayer),
    minBottomLayer_(minBottomLayer),
    maxBottomLayer_(maxBottomLayer),
    diameter_(diameter)
{}


Foam::engineValve::engineValve
(
    const word& name,
    const fvMeshMover& meshMover,
    const dictionary& dict
)
:
    name_(name),
    meshMover_(refCast<const fvMeshMovers::engine>(meshMover)),
    csPtr_
    (
        coordinateSystem::New
        (
            meshMover_.mesh(),
            dict.subDict("coordinateSystem")
        )
    ),
    bottomPatch_(dict.lookup("bottomPatch"), meshMover_.mesh().boundaryMesh()),
    poppetPatch_(dict.lookup("poppetPatch"), meshMover_.mesh().boundaryMesh()),
    stemPatch_(dict.lookup("stemPatch"), meshMover_.mesh().boundaryMesh()),
    curtainInPortPatch_
    (
        dict.lookup("curtainInPortPatch"),
        meshMover_.mesh().boundaryMesh()
    ),
    curtainInCylinderPatch_
    (
        dict.lookup("curtainInCylinderPatch"),
        meshMover_.mesh().boundaryMesh()
    ),
    detachInCylinderPatch_
    (
        dict.lookup("detachInCylinderPatch"),
        meshMover_.mesh().boundaryMesh()
    ),
    detachInPortPatch_
    (
        dict.lookup("detachInPortPatch"),
        meshMover_.mesh().boundaryMesh()
    ),
    detachFaces_(dict.lookup("detachFaces")),
    liftProfile_("liftProfile", dict),
    liftProfileStart_(min(liftProfile_.x())),
    liftProfileEnd_(max(liftProfile_.x())),
    minLift_(dict.lookup<scalar>("minLift")),
    minTopLayer_(dict.lookup<scalar>("minTopLayer")),
    maxTopLayer_(dict.lookup<scalar>("maxTopLayer")),
    minBottomLayer_(dict.lookup<scalar>("minBottomLayer")),
    maxBottomLayer_(dict.lookup<scalar>("maxBottomLayer")),
    diameter_(dict.lookup<scalar>("diameter"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::engineValve::lift(const scalar theta) const
{
    return liftProfile_.value(adjustCrankAngle(theta));
}


bool Foam::engineValve::isOpen() const
{
    return lift(meshMover_.theta()) >= minLift_;
}


Foam::scalar Foam::engineValve::curLift() const
{
    return max
    (
        lift(meshMover_.theta()),
        minLift_
    );
}


Foam::scalar Foam::engineValve::curVelocity() const
{
    return
       -(
             curLift()
           - max
             (
                 lift(meshMover_.theta() - meshMover_.deltaTheta()),
                 minLift_
             )
       )/(meshMover_.mesh().time().deltaTValue() + vSmall);
}


Foam::labelList Foam::engineValve::movingPatchIndices() const
{
    labelList mpIDs(2);
    label nMpIDs = 0;

    if (bottomPatch_.active())
    {
        mpIDs[nMpIDs] = bottomPatch_.index();
        nMpIDs++;
    }

    if (poppetPatch_.active())
    {
        mpIDs[nMpIDs] = poppetPatch_.index();
        nMpIDs++;
    }

    mpIDs.setSize(nMpIDs);

    return mpIDs;
}


void Foam::engineValve::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK;

    cs().writeDict(os);

    os  << "bottomPatch " << bottomPatch_.name() << token::END_STATEMENT << nl
        << "poppetPatch " << poppetPatch_.name() << token::END_STATEMENT << nl
        << "stemPatch " << stemPatch_.name() << token::END_STATEMENT << nl
        << "curtainInPortPatch " << curtainInPortPatch_.name()
        << token::END_STATEMENT << nl
        << "curtainInCylinderPatch " << curtainInCylinderPatch_.name()
        << token::END_STATEMENT << nl
        << "detachInCylinderPatch " << detachInCylinderPatch_.name()
        << token::END_STATEMENT << nl
        << "detachInPortPatch " << detachInPortPatch_.name()
        << token::END_STATEMENT << nl
        << "detachFaces " << detachFaces_ << token::END_STATEMENT << nl
        << "liftProfile " << nl << token::BEGIN_BLOCK
        << liftProfile_ << token::END_BLOCK << token::END_STATEMENT << nl
        << "minLift " << minLift_ << token::END_STATEMENT << nl
        << "minTopLayer " << minTopLayer_ << token::END_STATEMENT << nl
        << "maxTopLayer " << maxTopLayer_ << token::END_STATEMENT << nl
        << "minBottomLayer " << minBottomLayer_ << token::END_STATEMENT << nl
        << "maxBottomLayer " << maxBottomLayer_ << token::END_STATEMENT << nl
        << "diameter " << diameter_ << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //
