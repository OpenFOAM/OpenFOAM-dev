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

#include "displacementLaplacian_fvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "OFstream.H"
#include "meshTools.H"
#include "polyTopoChangeMap.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMotionSolvers
{
    defineTypeNameAndDebug(displacementLaplacian, 0);

    addToRunTimeSelectionTable
    (
        fvMeshMover,
        displacementLaplacian,
        fvMesh
    );

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        displacementLaplacian,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::displacementLaplacian::displacementLaplacian
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    fvMotionSolver(mesh),
    pointMeshMovers::displacement(mesh, dict, typeName),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMotionSolver::mesh(),
        dimensionedVector
        (
            "cellDisplacement",
            pointDisplacement_.dimensions(),
            Zero
        ),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
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
        Info<< "displacementLaplacian:" << nl
            << "    diffusivity       : " << diffusivityPtr_().type() << nl
            << "    frozenPoints zone : " << frozenPointsZone_ << endl;
    }

    if (io.headerOk())
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                io,
                pointMesh::New(fvMotionSolver::mesh())
            )
        );

        if (debug)
        {
            Info<< "displacementLaplacian :"
                << " Read pointVectorField "
                << io.name()
                << " to be used for boundary conditions on points."
                << nl
                << "Boundary conditions:"
                << pointLocation_().boundaryField().types() << endl;
        }
    }
}


Foam::fvMotionSolvers::displacementLaplacian::displacementLaplacian
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    displacementLaplacian(mesh.poly(), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::displacementLaplacian::~displacementLaplacian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::motionDiffusivity&
Foam::fvMotionSolvers::displacementLaplacian::diffusivity()
{
    if (!diffusivityPtr_.valid())
    {
        diffusivityType_.rewind();
        diffusivityPtr_ = motionDiffusivity::New
        (
            mesh(),
            diffusivityType_
        );
    }
    return diffusivityPtr_();
}


Foam::tmp<Foam::pointField>
Foam::fvMotionSolvers::displacementLaplacian::newPoints()
{
    // The points have moved so before interpolation update
    // the pointMeshMover accordingly
    movePoints(mesh().points());

    diffusivity().correct();
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    Foam::solve
    (
        fvm::laplacian
        (
            diffusivity().operator()(),
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
            Info<< "displacementLaplacian : applying "
                << " boundary conditions on " << pointLocation_().name()
                << " to new point location."
                << endl;
        }

        pointLocation_().primitiveFieldRef() =
            points0()
          + pointDisplacement_.primitiveField();

        pointLocation_().correctBoundaryConditions();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = mesh().pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                pointLocation_()[pz[i]] = points0()[pz[i]];
            }
        }

        twoDCorrectPoints(pointLocation_().primitiveFieldRef());

        return tmp<pointField>(pointLocation_().primitiveField());
    }
    else
    {
        tmp<pointField> tcurPoints
        (
            points0() + pointDisplacement_.primitiveField()
        );
        pointField& curPoints = tcurPoints.ref();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = mesh().pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                curPoints[pz[i]] = points0()[pz[i]];
            }
        }

        twoDCorrectPoints(curPoints);

        return tcurPoints;
    }
}


void Foam::fvMotionSolvers::displacementLaplacian::topoChange
(
    const polyTopoChangeMap& map
)
{
    pointMeshMovers::displacement::topoChange(map);
    diffusivityPtr_.clear();
}


void Foam::fvMotionSolvers::displacementLaplacian::mapMesh
(
    const polyMeshMap& map
)
{
    pointMeshMovers::displacement::mapMesh(map);
    diffusivityPtr_.clear();
}


// ************************************************************************* //
