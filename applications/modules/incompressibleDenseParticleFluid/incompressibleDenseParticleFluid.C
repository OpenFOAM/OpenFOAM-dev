/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "incompressibleDenseParticleFluid.H"
#include "surfaceInterpolate.H"
#include "fvMatrix.H"
#include "fvcFlux.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleDenseParticleFluid, 0);
    addToRunTimeSelectionTable
    (
        solver,
        incompressibleDenseParticleFluid,
        fvMesh
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::incompressibleDenseParticleFluid::correctCoNum()
{
    fluidSolver::correctCoNum(phic);
}


void Foam::solvers::incompressibleDenseParticleFluid::continuityErrors()
{
    fluidSolver::continuityErrors(phic);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleDenseParticleFluid::
incompressibleDenseParticleFluid
(
    fvMesh& mesh
)
:
    fluidSolver(mesh),

    continuousPhaseName
    (
        IOdictionary
        (
            IOobject
            (
                "physicalProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ
            )
        ).lookup("continuousPhaseName")
    ),

    p_
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    pressureReference(p_, pimple.dict()),

    g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    Uc_
    (
        IOobject
        (
            IOobject::groupName("U", continuousPhaseName),
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    phic_
    (
        IOobject
        (
            IOobject::groupName("phi", continuousPhaseName),
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(Uc_) & mesh.Sf()
    ),

    viscosity(viscosityModel::New(mesh)),

    rhoc
    (
        IOobject
        (
            IOobject::groupName("rho", continuousPhaseName),
            runTime.name(),
            mesh
        ),
        mesh,
        dimensionedScalar
        (
            IOobject::groupName("rho", continuousPhaseName),
            dimDensity,
            viscosity->lookup
            (
                IOobject::groupName("rho", continuousPhaseName)
            )
        )
    ),

    muc
    (
        IOobject
        (
            IOobject::groupName("mu", continuousPhaseName),
            runTime.name(),
            mesh
        ),
        rhoc*viscosity->nu()
    ),

    alphac_
    (
        IOobject
        (
            IOobject::groupName("alpha", continuousPhaseName),
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    ),

    alphacMin
    (
        1 - mesh.solution().solverDict(alphac_.name()).lookup<scalar>("max")
    ),

    alphacf("alphacf", fvc::interpolate(alphac_)),

    alphaPhic
    (
        IOobject::groupName("alphaPhi", continuousPhaseName),
        alphacf*phic_
    ),

    momentumTransport
    (
        phaseIncompressible::momentumTransportModel::New
        (
            alphac_,
            Uc_,
            alphaPhic,
            phic_,
            viscosity
        )
    ),

    clouds
    (
        parcelClouds::New(mesh, rhoc, Uc_, muc, g)
    ),

    p(p_),
    Uc(Uc_),
    phic(phic_),
    alphac(alphac_)
{
    mesh.schemes().setFluxRequired(p.name());

    momentumTransport->validate();

    // Update alphac from the particle locations
    alphac_ = max(1 - clouds.alpha(), alphacMin);
    alphac_.correctBoundaryConditions();
    alphacf = fvc::interpolate(alphac);
    alphaPhic = alphacf*phic;

    correctCoNum();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleDenseParticleFluid::
~incompressibleDenseParticleFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleDenseParticleFluid::preSolve()
{
    if (mesh.dynamic() && !Ucf.valid())
    {
        Info<< "Constructing face momentum Ucf" << endl;

        // Ensure the U BCs are up-to-date before constructing Ucf
        Uc_.correctBoundaryConditions();

        Ucf = new surfaceVectorField
        (
            IOobject
            (
                "Ucf",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(Uc)
        );
    }

    // Store the particle positions
    if (mesh.topoChanging() || mesh.distributing())
    {
        clouds.storeGlobalPositions();
    }

    fvModels().preUpdateMesh();

    correctCoNum();

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}


void Foam::solvers::incompressibleDenseParticleFluid::prePredictor()
{
    if (pimple.firstIter())
    {
        clouds.evolve();

        // Update continuous phase volume fraction field
        alphac_ = max(1 - clouds.alpha(), alphacMin);
        alphac_.correctBoundaryConditions();
        alphacf = fvc::interpolate(alphac);

        // ... and continuous phase volumetric flux
        alphaPhic = alphacf*phic;

        // Update the continuous phase dynamic viscosity
        muc = rhoc*viscosity->nu();

        Fd = new volVectorField
        (
            IOobject
            (
                "Fd",
                runTime.name(),
                mesh
            ),
            mesh,
            dimensionedVector(dimAcceleration, Zero),
            zeroGradientFvPatchVectorField::typeName
        );

        Dc = new volScalarField
        (
            IOobject
            (
                "Dc",
                runTime.name(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero),
            zeroGradientFvPatchVectorField::typeName
        );

        const fvVectorMatrix cloudSU(clouds.SU(Uc));

        Fd().primitiveFieldRef() = -cloudSU.source()/mesh.V()/rhoc;
        Fd().correctBoundaryConditions();

        Dc().primitiveFieldRef() = -cloudSU.diag()/mesh.V()/rhoc;
        Dc().correctBoundaryConditions();

        Dcf = fvc::interpolate(Dc()).ptr();

        phid =
        (
            fvc::flux(Fd())
           /(Dcf() + dimensionedScalar(Dc().dimensions(), small))
        ).ptr();
    }

    if (pimple.predictTransport())
    {
        momentumTransport->predict();
    }
}


void Foam::solvers::incompressibleDenseParticleFluid::thermophysicalPredictor()
{}


void Foam::solvers::incompressibleDenseParticleFluid::pressureCorrector()
{
    while (pimple.correct())
    {
        correctPressure();
    }

    tUcEqn.clear();
}


void Foam::solvers::incompressibleDenseParticleFluid::postCorrector()
{
    if (pimple.correctTransport())
    {
        viscosity->correct();
        momentumTransport->correct();
    }
}


void Foam::solvers::incompressibleDenseParticleFluid::postSolve()
{
    Fd.clear();
    Dc.clear();
    Dcf.clear();
    phid.clear();
}


// ************************************************************************* //
