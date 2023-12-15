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

#include "displacementComponentLaplacianFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "polyTopoChangeMap.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementComponentLaplacianFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementComponentLaplacianFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementComponentLaplacianFvMotionSolver::
displacementComponentLaplacianFvMotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    componentDisplacementMotionSolver(name, mesh, dict, type()),
    fvMotionSolver(mesh),
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
        fvMesh_,
        dimensionedScalar(pointDisplacement_.dimensions(), 0),
        cellMotionBoundaryTypes<scalar>(pointDisplacement_.boundaryField())
    ),
    pointLocation_(nullptr),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findIndex(coeffDict().lookup("frozenPointsZone"))
      : -1
    )
{
    typeIOobject<pointVectorField> io
    (
        "pointLocation",
        fvMesh_.time().name(),
        fvMesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (debug)
    {
        Info<< "displacementComponentLaplacianFvMotionSolver:" << nl
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
                    fvMesh_.time().name(),
                    fvMesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                pointMesh::New(fvMesh_)
            )
        );

        if (debug)
        {
            Info<< "displacementComponentLaplacianFvMotionSolver :"
                << " Read pointVectorField "
                << pointLocation_().name()
                << " to be used for boundary conditions on points."
                << nl
                << "Boundary conditions:"
                << pointLocation_().boundaryField().types() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementComponentLaplacianFvMotionSolver::
~displacementComponentLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::displacementComponentLaplacianFvMotionSolver::curPoints() const
{
    volPointInterpolation::New(fvMesh_).interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    if (pointLocation_.valid())
    {
        if (debug)
        {
            Info<< "displacementComponentLaplacianFvMotionSolver : applying "
                << " boundary conditions on " << pointLocation_().name()
                << " to new point location."
                << endl;
        }

        // Apply pointLocation_ b.c. to mesh points.

        pointLocation_().primitiveFieldRef() = fvMesh_.points();

        pointLocation_().primitiveFieldRef().replace
        (
            cmpt_,
            points0_ + pointDisplacement_.primitiveField()
        );

        pointLocation_().correctBoundaryConditions();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

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
        tmp<pointField> tcurPoints(new pointField(fvMesh_.points()));
        pointField& curPoints = tcurPoints.ref();

        curPoints.replace
        (
            cmpt_,
            points0_ + pointDisplacement_.primitiveField()
        );

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

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


void Foam::displacementComponentLaplacianFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the motionSolver accordingly
    movePoints(fvMesh_.points());

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
}


void Foam::displacementComponentLaplacianFvMotionSolver::topoChange
(
    const polyTopoChangeMap& map
)
{
    componentDisplacementMotionSolver::topoChange(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityPtr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity")
    );
}


void Foam::displacementComponentLaplacianFvMotionSolver::mapMesh
(
    const polyMeshMap& map
)
{
    componentDisplacementMotionSolver::mapMesh(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityPtr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity")
    );
}


// ************************************************************************* //
