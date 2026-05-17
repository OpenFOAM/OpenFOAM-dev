/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "velocityComponentLaplacian_fvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMotionSolvers
{
    defineTypeNameAndDebug(velocityComponentLaplacian, 0);

    addToRunTimeSelectionTable
    (
        fvMeshMover,
        velocityComponentLaplacian,
        fvMesh
    );

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        velocityComponentLaplacian,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::velocityComponentLaplacian::velocityComponentLaplacian
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    fvMotionSolver(mesh),
    pointMeshMovers::velocityComponent(mesh, dict, typeName),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU" + cmptName_,
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMotionSolver::mesh(),
        dimensionedScalar(pointMotionU_.dimensions(), 0),
        cellMotionBoundaryTypes<scalar>(pointMotionU_.boundaryField())
    ),
    diffusivityType_(dict.lookup("diffusivity")),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMotionSolver::mesh(), diffusivityType_)
    )
{}


Foam::fvMotionSolvers::velocityComponentLaplacian::velocityComponentLaplacian
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    velocityComponentLaplacian(mesh.poly(), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::velocityComponentLaplacian::~velocityComponentLaplacian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::fvMotionSolvers::velocityComponentLaplacian::newPoints()
{
    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(mesh().points());

    diffusivityPtr_->correct();
    pointMotionU_.boundaryFieldRef().updateCoeffs();

    Foam::solve
    (
        fvm::laplacian
        (
            diffusivityPtr_->operator()(),
            cellMotionU_,
            "laplacian(diffusivity,cellMotionU)"
        )
    );

    volPointInterpolation::New(mesh()).interpolate
    (
        cellMotionU_,
        pointMotionU_
    );

    tmp<pointField> tcurPoints(new pointField(mesh().points()));

    tcurPoints.ref().replace
    (
        cmpt_,
        tcurPoints().component(cmpt_)
      + mesh().time().deltaTValue()
       *pointMotionU_.primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::fvMotionSolvers::velocityComponentLaplacian::topoChange
(
    const polyTopoChangeMap& map
)
{
    pointMeshMovers::velocityComponent::topoChange(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityType_.rewind();
    diffusivityPtr_ = motionDiffusivity::New
    (
        mesh(),
        diffusivityType_
    );
}


void Foam::fvMotionSolvers::velocityComponentLaplacian::mapMesh
(
    const polyMeshMap& map
)
{
    pointMeshMovers::velocityComponent::mapMesh(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityType_.rewind();
    diffusivityPtr_ = motionDiffusivity::New
    (
        mesh(),
        diffusivityType_
    );
}


// ************************************************************************* //
