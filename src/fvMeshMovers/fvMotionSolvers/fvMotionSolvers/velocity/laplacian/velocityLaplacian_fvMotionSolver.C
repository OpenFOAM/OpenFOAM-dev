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

#include "velocityLaplacian_fvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMotionSolvers
{
    defineTypeNameAndDebug(velocityLaplacian, 0);

    addToRunTimeSelectionTable
    (
        fvMeshMover,
        velocityLaplacian,
        fvMesh
    );

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        velocityLaplacian,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::velocityLaplacian::velocityLaplacian
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    fvMotionSolver(mesh),
    pointMeshMovers::velocity(mesh, dict, typeName),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMotionSolver::mesh(),
        dimensionedVector
        (
            "cellMotionU",
            pointMotionU_.dimensions(),
            Zero
        ),
        cellMotionBoundaryTypes<vector>(pointMotionU_.boundaryField())
    ),
    diffusivityType_(dict.lookup("diffusivity")),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMotionSolver::mesh(), diffusivityType_)
    )
{}


Foam::fvMotionSolvers::velocityLaplacian::velocityLaplacian
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    velocityLaplacian(mesh.poly(), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::velocityLaplacian::~velocityLaplacian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::fvMotionSolvers::velocityLaplacian::newPoints()
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

    tmp<pointField> tcurPoints
    (
        mesh().points()
      + mesh().time().deltaTValue()*pointMotionU_.primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


//void Foam::fvMotionSolvers::velocityLaplacian::movePoints(const pointField& p)
//{
//    // Movement of pointMesh and volPointInterpolation already
//    // done by polyMesh,fvMesh
//}


void Foam::fvMotionSolvers::velocityLaplacian::topoChange
(
    const polyTopoChangeMap& map
)
{
    pointMeshMovers::velocity::topoChange(map);

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


void Foam::fvMotionSolvers::velocityLaplacian::mapMesh
(
    const polyMeshMap& map
)
{
    pointMeshMovers::velocity::mapMesh(map);

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
