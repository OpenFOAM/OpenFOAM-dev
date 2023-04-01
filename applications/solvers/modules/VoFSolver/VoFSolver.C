/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "VoFSolver.H"
#include "localEulerDdtScheme.H"
#include "linear.H"
#include "fvcDiv.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(VoFSolver, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::VoFSolver::continuityErrors()
{
    fluidSolver::continuityErrors(phi);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::solvers::VoFSolver::setrAU(const fvVectorMatrix& UEqn)
{
    if (rAU.valid())
    {
        rAU() = 1.0/UEqn.A();
    }
    else
    {
        rAU = (1.0/UEqn.A()).ptr();
    }
}


void Foam::solvers::VoFSolver::clearrAU()
{
    if (!(correctPhi || mesh.topoChanging()))
    {
        rAU.clear();
    }
}


void Foam::solvers::VoFSolver::correctCoNum()
{
    fluidSolver::correctCoNum(phi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::VoFSolver::VoFSolver
(
    fvMesh& mesh,
    autoPtr<VoFMixture> mixturePtr
)
:
    fluidSolver(mesh),

    mixture_(mixturePtr),
    mixture(mixture_()),

    divAlphaName("div(phi,alpha)"),

    U_
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    phi
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U_) & mesh.Sf()
    ),

    buoyancy(mesh),

    p_rgh(buoyancy.p_rgh),

    rho(mixture.rho()),

    rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    ),

    MRF(mesh),

    U(U_)
{
    mesh.schemes().setFluxRequired(p_rgh.name());

    if (mesh.dynamic() || MRF.size())
    {
        Info<< "Constructing face momentum Uf" << endl;

        // Ensure the U BCs are up-to-date before constructing Uf
        U_.correctBoundaryConditions();

        Uf = new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U_)
        );
    }

    if (LTS)
    {
        Info<< "Using LTS" << endl;

        trDeltaT = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime.name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::VoFSolver::~VoFSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::VoFSolver::maxDeltaT() const
{
    const scalar maxAlphaCo =
        runTime.controlDict().lookup<scalar>("maxAlphaCo");

    scalar deltaT = fluidSolver::maxDeltaT();

    if (alphaCoNum > small)
    {
        deltaT = min(deltaT, maxAlphaCo/alphaCoNum*runTime.deltaTValue());
    }

    return deltaT;
}


void Foam::solvers::VoFSolver::preSolve()
{
    // Read the controls
    readControls();

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        setRDeltaT();
    }

    fvModels().preUpdateMesh();

    // Store divU from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    if (correctPhi && divergent())
    {
        // Construct and register divU for mapping
        divU = new volScalarField
        (
            "divU0",
            fvc::div(fvc::absolute(phi, U))
        );
    }

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}


void Foam::solvers::VoFSolver::prePredictor()
{
    fvModels().correct();
}


void Foam::solvers::VoFSolver::postSolve()
{
    divU.clear();
}


// ************************************************************************* //
