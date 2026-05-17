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

#include "displacementComponentLaplacian_fvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "polyTopoChangeMap.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMotionSolvers
{
    defineTypeNameAndDebug(displacementComponentLaplacian, 0);

    addToRunTimeSelectionTable
    (
        fvMeshMover,
        displacementComponentLaplacian,
        fvMesh
    );

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        displacementComponentLaplacian,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::displacementComponentLaplacian::
displacementComponentLaplacian
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    fvMotionSolver(mesh),
    pointMeshMovers::displacementComponent(mesh, dict, type()),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement" + cmptName_,
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMotionSolver::mesh(),
        dimensionedScalar(pointDisplacement_.dimensions(), 0),
        cellMotionBoundaryTypes<scalar>(pointDisplacement_.boundaryField())
    ),
    pointLocation_(nullptr),
    diffusivityType_(dict.lookup("diffusivity")),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMotionSolver::mesh(), diffusivityType_)
    ),
    frozenPointsZone_
    (
        dict.found("frozenPointsZone")
      ? fvMotionSolver::mesh().pointZones().findIndex
        (
            dict.lookup("frozenPointsZone")
        )
      : -1
    )
{
    typeIOobject<pointVectorField> io
    (
        "pointLocation",
        fvMotionSolver::mesh().time().name(),
        fvMotionSolver::mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (debug)
    {
        Info<< "displacementComponentLaplacian:" << nl
            << "    diffusivity       : " << diffusivityPtr_().type() << nl
            << "    frozenPoints zone : " << frozenPointsZone_ << endl;
    }

    if (io.headerOk())
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                IOobject
                (
                    "pointLocation",
                    fvMotionSolver::mesh().time().name(),
                    fvMotionSolver::mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                pointMesh::New(fvMotionSolver::mesh())
            )
        );

        if (debug)
        {
            Info<< "displacementComponentLaplacian :"
                << " Read pointVectorField "
                << pointLocation_().name()
                << " to be used for boundary conditions on points."
                << nl
                << "Boundary conditions:"
                << pointLocation_().boundaryField().types() << endl;
        }
    }
}


Foam::fvMotionSolvers::displacementComponentLaplacian::
displacementComponentLaplacian
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    displacementComponentLaplacian(mesh.poly(), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::displacementComponentLaplacian::
~displacementComponentLaplacian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::fvMotionSolvers::displacementComponentLaplacian::newPoints()
{
    // The points have moved so before interpolation update
    // the pointMeshMover accordingly
    movePoints(mesh().points());

    diffusivityPtr_->correct();
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    Foam::solve
    (
        fvm::laplacian
        (
            diffusivityPtr_->operator()(),
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )
    );

    volPointInterpolation::New(mesh()).interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    if (pointLocation_.valid())
    {
        if (debug)
        {
            Info<< "displacementComponentLaplacian : applying "
                << " boundary conditions on " << pointLocation_().name()
                << " to new point location."
                << endl;
        }

        // Apply pointLocation_ b.c. to mesh points.

        pointLocation_().primitiveFieldRef() = mesh().points();

        pointLocation_().primitiveFieldRef().replace
        (
            cmpt_,
            points0_ + pointDisplacement_.primitiveField()
        );

        pointLocation_().correctBoundaryConditions();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = mesh().pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                label pointi = pz[i];

                pointLocation_()[pointi][cmpt_] = points0_[pointi];
            }
        }

        twoDCorrectPoints(pointLocation_().primitiveFieldRef());

        return tmp<pointField>(pointLocation_().primitiveField());
    }
    else
    {
        tmp<pointField> tcurPoints(new pointField(mesh().points()));
        pointField& curPoints = tcurPoints.ref();

        curPoints.replace
        (
            cmpt_,
            points0_ + pointDisplacement_.primitiveField()
        );

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = mesh().pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                label pointi = pz[i];

                curPoints[pointi][cmpt_] = points0_[pointi];
            }
        }

        twoDCorrectPoints(curPoints);

        return tcurPoints;
    }
}


void Foam::fvMotionSolvers::displacementComponentLaplacian::topoChange
(
    const polyTopoChangeMap& map
)
{
    pointMeshMovers::displacementComponent::topoChange(map);

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


void Foam::fvMotionSolvers::displacementComponentLaplacian::mapMesh
(
    const polyMeshMap& map
)
{
    pointMeshMovers::displacementComponent::mapMesh(map);

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
