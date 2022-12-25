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

#include "VoFSolver.H"
#include "localEulerDdtScheme.H"
#include "CorrectPhi.H"
#include "geometricZeroField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(VoFSolver, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::VoFSolver::correctCoNum()
{
    fluidSolver::correctCoNum(phi);

    const scalarField sumPhi
    (
        interface.nearInterface()().primitiveField()
       *fvc::surfaceSum(mag(phi))().primitiveField()
    );

    alphaCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    const scalar meanAlphaCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    Info<< "Interface Courant Number mean: " << meanAlphaCoNum
        << " max: " << alphaCoNum << endl;
}


void Foam::solvers::VoFSolver::continuityErrors()
{
    fluidSolver::continuityErrors(phi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::VoFSolver::VoFSolver
(
    fvMesh& mesh,
    autoPtr<twoPhaseMixture> mixturePtr
)
:
    fluidSolver(mesh),

    mixture_(mixturePtr),
    mixture(mixture_()),

    alpha1(mixture.alpha1()),
    alpha2(mixture.alpha2()),

    alphaRestart
    (
        typeIOobject<surfaceScalarField>
        (
            IOobject::groupName("alphaPhi", alpha1.group()),
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ).headerOk()
    ),

    divAlphaName("div(phi,alpha)"),

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

    interface(mixture, alpha1, alpha2, U),

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

    alphaPhi1
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", alpha1.group()),
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        phi*fvc::interpolate(alpha1)
    ),

    MRF(mesh)
{
    // Read the controls
    read();

    mesh.schemes().setFluxRequired(p_rgh.name());
    mesh.schemes().setFluxRequired(alpha1.name());

    if (alphaRestart)
    {
        Info << "Restarting alpha" << endl;
    }

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

Foam::solvers::VoFSolver::~VoFSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::VoFSolver::maxDeltaT() const
{
    const scalar maxAlphaCo
    (
        runTime.controlDict().lookup<scalar>("maxAlphaCo")
    );

    const scalar deltaT = fluidSolver::maxDeltaT();

    if (alphaCoNum > small)
    {
        return min
        (
            deltaT,
            maxAlphaCo/(alphaCoNum + small)*runTime.deltaTValue()
        );
    }
    else
    {
        return deltaT;
    }
}


void Foam::solvers::VoFSolver::preSolve()
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

    fvModels().preUpdateMesh();

    // Store divU from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    if
    (
        correctPhi
     && divergent()
     && mesh.topoChanged())
    {
        // Construct and register divU for mapping
        divU = new volScalarField
        (
            "divU0",
            fvc::div(fvc::absolute(phi, U))
        );
    }

    // Update the mesh for topology change, mesh to mesh mapping
    const bool topoChanged = mesh.update();

    // Do not apply previous time-step mesh compression flux
    // if the mesh topology changed
    if (topoChanged)
    {
        talphaPhi1Corr0.clear();
    }
}


void Foam::solvers::VoFSolver::prePredictor()
{
    fvModels().correct();
    alphaPredictor();
}


void Foam::solvers::VoFSolver::postCorrector()
{}


void Foam::solvers::VoFSolver::postSolve()
{
    divU.clear();
}


// ************************************************************************* //
