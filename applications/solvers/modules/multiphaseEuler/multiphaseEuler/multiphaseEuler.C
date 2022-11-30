/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "multiphaseEuler.H"
#include "localEulerDdtScheme.H"
#include "surfaceFields.H"
#include "fvcDiv.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(multiphaseEuler, 0);
    addToRunTimeSelectionTable(solver, multiphaseEuler, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::multiphaseEuler::read()
{
    fluidSolver::read();

    faceMomentum =
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false);

    partialElimination =
        pimple.dict().lookupOrDefault<Switch>("partialElimination", false);

    nEnergyCorrectors =
        pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1);
}


void Foam::solvers::multiphaseEuler::correctCoNum()
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    forAll(phases, phasei)
    {
        sumPhi = max
        (
            sumPhi,
            fvc::surfaceSum(mag(phases[phasei].phi()))().primitiveField()
        );
    }

    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    const scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::multiphaseEuler::multiphaseEuler(fvMesh& mesh)
:
    fluidSolver(mesh),

    faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    ),

    partialElimination
    (
        pimple.dict().lookupOrDefault<Switch>("partialElimination", false)
    ),

    nEnergyCorrectors
    (
        pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1)
    ),

    trDeltaT
    (
        LTS
      ? new volScalarField
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
      : nullptr
    ),

    trDeltaTf
    (
        LTS && faceMomentum
      ? new surfaceScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTfName,
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, 1)
        )
      : nullptr
    ),

    buoyancy(mesh),

    fluidPtr(phaseSystem::New(mesh)),

    fluid(fluidPtr()),

    phases(fluid.phases()),

    phi(fluid.phi()),

    p(phases[0].thermoRef().p()),

    p_rgh(buoyancy.p_rgh),

    pressureReference
    (
        p,
        p_rgh,
        pimple.dict(),
        fluid.incompressible()
    ),

    MRF(fluid.MRF())
{
    // Read the controls
    read();

    mesh.schemes().setFluxRequired(p_rgh.name());

    if (transient())
    {
        correctCoNum();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::multiphaseEuler::~multiphaseEuler()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multiphaseEuler::preSolve()
{
    // Read the controls
    read();

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        setRDeltaT();
    }

    // Store divU from the previous mesh so that it can be
    // mapped and used in correctPhi to ensure the corrected phi
    // has the same divergence
    if (correctPhi)
    {
        // Construct and register divU for mapping
        divU = new volScalarField
        (
            "divU0",
            fvc::div
            (
                fvc::absolute(phi, fluid.movingPhases()[0].U())
            )
        );
    }

    fvModels().preUpdateMesh();

    // Update the mesh for topology change, mesh to mesh mapping
    mesh.update();
}


void Foam::solvers::multiphaseEuler::prePredictor()
{
    if (pimple.models())
    {
        fvModels().correct();
    }

    if (pimple.thermophysics() || pimple.flow())
    {
        fluid.solve(rAUs, rAUfs);
        fluid.correct();
        fluid.correctContinuityError();
    }

    if (pimple.thermophysics())
    {
        compositionPredictor();
    }
}


void Foam::solvers::multiphaseEuler::momentumTransportCorrector()
{
    if (pimple.flow() && pimple.transportCorr())
    {
        fluid.correctTurbulence();
    }
}


void Foam::solvers::multiphaseEuler::thermophysicalTransportCorrector()
{}


void Foam::solvers::multiphaseEuler::postSolve()
{
    divU.clear();
}


// ************************************************************************* //
