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

#include "incompressibleFluid.H"
#include "localEulerDdtScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleFluid, 0);
    addToRunTimeSelectionTable(solver, incompressibleFluid, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::incompressibleFluid::correctCoNum()
{
    fluidSolver::correctCoNum(phi);
}


void Foam::solvers::incompressibleFluid::continuityErrors()
{
    fluidSolver::continuityErrors(phi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleFluid::incompressibleFluid(fvMesh& mesh)
:
    fluidSolver(mesh),

    p
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

    pressureReference(p, pimple.dict()),

    U
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
        linearInterpolate(U) & mesh.Sf()
    ),

    viscosity(viscosityModel::New(mesh)),

    momentumTransport
    (
        incompressible::momentumTransportModel::New
        (
            U,
            phi,
            viscosity
        )
    ),

    MRF(mesh)
{
    mesh.schemes().setFluxRequired(p.name());

    momentumTransport->validate();

    if (mesh.dynamic())
    {
        Info<< "Constructing face momentum Uf" << endl;

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
            fvc::interpolate(U)
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

Foam::solvers::incompressibleFluid::~incompressibleFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleFluid::preSolve()
{
    // Read the controls
    read();

    fvModels().preUpdateMesh();

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        setRDeltaT();
    }

    // Update the mesh for topology change, mesh to mesh mapping
    mesh.update();
}


void Foam::solvers::incompressibleFluid::prePredictor()
{
    fvModels().correct();
}


void Foam::solvers::incompressibleFluid::thermophysicalPredictor()
{}


void Foam::solvers::incompressibleFluid::pressureCorrector()
{
    while (pimple.correct())
    {
        correctPressure();
    }

    tUEqn.clear();
}


void Foam::solvers::incompressibleFluid::momentumTransportCorrector()
{
    if (pimple.transportCorr())
    {
        viscosity->correct();
        momentumTransport->correct();
    }
}


void Foam::solvers::incompressibleFluid::thermophysicalTransportCorrector()
{}


void Foam::solvers::incompressibleFluid::postSolve()
{}


// ************************************************************************* //
