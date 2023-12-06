/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "incompressibleDriftFlux.H"
#include "fvCorrectPhi.H"
#include "addToRunTimeSelectionTable.H"

#include "fvmDdt.H"

#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleDriftFlux, 0);
    addToRunTimeSelectionTable(solver, incompressibleDriftFlux, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::incompressibleDriftFlux::correctCoNum()
{
    VoFSolver::correctCoNum();
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::solvers::incompressibleDriftFlux::setInterfaceRDeltaT
(
    volScalarField& rDeltaT
)
{}


void Foam::solvers::incompressibleDriftFlux::correctInterface()
{}


Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::incompressibleDriftFlux::surfaceTensionForce() const
{
    return surfaceScalarField::New
    (
        "surfaceTensionForce",
        mesh,
        dimensionedScalar(dimForce/dimVolume, 0)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleDriftFlux::incompressibleDriftFlux(fvMesh& mesh)
:
    twoPhaseSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new incompressibleDriftFluxMixture(mesh))
    ),

    mixture
    (
        refCast<incompressibleDriftFluxMixture>(twoPhaseSolver::mixture)
       .initialise(U)
    ),

    p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*buoyancy.gh
    ),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict()
    ),

    relativeVelocity
    (
        relativeVelocityModel::New(mixture, mixture, buoyancy.g)
    ),

    momentumTransport
    (
        compressible::momentumTransportModel::New(rho, U, rhoPhi, mixture)
    )
{
    if (transient())
    {
        correctCoNum();
    }

    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }

    if (!runTime.restart() || !divergent())
    {
        correctUphiBCs(U_, phi_, true);

        fv::correctPhi
        (
            phi_,
            U,
            p_rgh,
            rAU,
            autoPtr<volScalarField>(),
            pressureReference(),
            pimple
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleDriftFlux::~incompressibleDriftFlux()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::incompressibleDriftFlux::maxDeltaT() const
{
    return fluidSolver::maxDeltaT();
}


void Foam::solvers::incompressibleDriftFlux::prePredictor()
{
    VoFSolver::prePredictor();
    alphaPredictor();

    // Apply the diffusion term separately to allow implicit solution
    // and boundedness of the explicit advection
    {
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(alpha1) - fvc::ddt(alpha1)
          - fvm::laplacian(momentumTransport->nut(), alpha1)
        );

        alpha1Eqn.solve(alpha1.name() + "Diffusion");

        alphaPhi1 += alpha1Eqn.flux();

        alpha2 = scalar(1) - alpha1;
        alphaPhi2 = phi - alphaPhi1;

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }

    mixture.correct();

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    // Calculate the mass-flux
    rhoPhi = alphaPhi1*rho1 + alphaPhi2*rho2;

    relativeVelocity->correct();

    if (pimple.predictTransport())
    {
        momentumTransport->predict();
    }
}


void Foam::solvers::incompressibleDriftFlux::pressureCorrector()
{
    incompressiblePressureCorrector(p);
}


void Foam::solvers::incompressibleDriftFlux::thermophysicalPredictor()
{}


void Foam::solvers::incompressibleDriftFlux::postCorrector()
{
    if (pimple.correctTransport())
    {
        momentumTransport->correct();
    }
}


// ************************************************************************* //
