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

#include "isothermalFluid.H"
#include "localEulerDdtScheme.H"
#include "hydrostaticInitialisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(isothermalFluid, 0);
    addToRunTimeSelectionTable(solver, isothermalFluid, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::isothermalFluid::correctCoNum()
{
    fluidSolver::correctCoNum(rho, phi);
}


void Foam::solvers::isothermalFluid::continuityErrors()
{
    fluidSolver::continuityErrors(rho, thermo.rho(), phi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::isothermalFluid::isothermalFluid
(
    fvMesh& mesh,
    autoPtr<fluidThermo> thermoPtr
)
:
    fluidSolver(mesh),

    thermo_(thermoPtr),
    thermo(thermo_()),

    p(thermo.p()),

    rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    ),

    dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(p.dimensions()/dimTime, 0)
    ),

    buoyancy(buoyancy::New(mesh)),

    p_rgh(buoyancy.valid() ? buoyancy->p_rgh : p),

    pressureReference
    (
        p,
        p_rgh,
        pimple.dict(),
        thermo.incompressible()
    ),

    U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
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
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho*U) & mesh.Sf()
    ),

    K("K", 0.5*magSqr(U)),

    turbulence
    (
        compressible::momentumTransportModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    ),

    initialMass(fvc::domainIntegrate(rho)),

    MRF(mesh)
{
    // Read the controls
    read();

    thermo.validate("isothermalFluid", "h", "e");
    mesh.schemes().setFluxRequired(p.name());
    turbulence->validate();

    if (buoyancy.valid())
    {
        hydrostaticInitialisation
        (
            p_rgh,
            p,
            rho,
            U,
            buoyancy->gh,
            buoyancy->ghf,
            buoyancy->pRef,
            thermo,
            pimple.dict()
        );
    }

    if (mesh.dynamic())
    {
        Info<< "Constructing face momentum rhoUf" << endl;

        rhoUf = new surfaceVectorField
        (
            IOobject
            (
                "rhoUf",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(rho*U)
        );
    }

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        Info<< "Using LTS" << endl;

        trDeltaT = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime.timeName(),
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


Foam::solvers::isothermalFluid::isothermalFluid(fvMesh& mesh)
:
    isothermalFluid(mesh, fluidThermo::New(mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::isothermalFluid::~isothermalFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::isothermalFluid::preSolve()
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

    // Store divrhoU from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    if (correctPhi)
    {
        divrhoU = new volScalarField
        (
            "divrhoU",
            fvc::div(fvc::absolute(phi, rho, U))
        );
    }

    fvModels().preUpdateMesh();

    // Store momentum to set rhoUf for introduced faces
    if (rhoUf.valid())
    {
        rhoU = new volVectorField("rhoU", rho*U);
    }

    // Update the mesh for topology change, mesh to mesh mapping
    mesh.update();
}


void Foam::solvers::isothermalFluid::thermophysicalPredictor()
{
    thermo.correct();
}


void Foam::solvers::isothermalFluid::pressureCorrector()
{
    while (pimple.correct())
    {
        if (buoyancy.valid())
        {
            correctBuoyantPressure();
        }
        else
        {
            correctPressure();
        }
    }

    tUEqn.clear();
}


void Foam::solvers::isothermalFluid::momentumTransportCorrector()
{
    if (pimple.transportCorr())
    {
        turbulence->correct();
    }
}


void Foam::solvers::isothermalFluid::thermophysicalTransportCorrector()
{}


void Foam::solvers::isothermalFluid::postSolve()
{
    rhoU.clear();
    divrhoU.clear();

    if (!mesh.schemes().steady())
    {
        rho = thermo.rho();
    }
}


// ************************************************************************* //
