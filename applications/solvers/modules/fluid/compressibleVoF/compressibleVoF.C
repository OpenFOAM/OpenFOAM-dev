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

#include "compressibleVoF.H"
#include "localEulerDdtScheme.H"
#include "hydrostaticInitialisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(compressibleVoF, 0);
    addToRunTimeSelectionTable(solver, compressibleVoF, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::compressibleVoF::correctCoNum()
{
    fluidSolver::correctCoNum(phi);

    alphaCoNum = 0;
    scalar meanAlphaCoNum = 0;

    if (mesh.nInternalFaces())
    {
        const scalarField sumPhi
        (
            mixture.nearInterface()().primitiveField()
            *fvc::surfaceSum(mag(phi))().primitiveField()
        );

        alphaCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

        meanAlphaCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }

    Info<< "Interface Courant Number mean: " << meanAlphaCoNum
        << " max: " << alphaCoNum << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::compressibleVoF::compressibleVoF(fvMesh& mesh)
:
    fluidSolver(mesh),

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
        linearInterpolate(U) & mesh.Sf()
    ),

    mixture(U, phi),

    alpha1(mixture.alpha1()),

    alphaRestart
    (
        typeIOobject<surfaceScalarField>
        (
            IOobject::groupName("alphaPhi", alpha1.group()),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ).headerOk()
    ),

    buoyancy(mesh),

    p_rgh(buoyancy.p_rgh),

    rho(mixture.rho()),

    dgdt
    (
        IOobject
        (
            "dgdt",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*fvc::div(phi)
    ),

    pressureReference
    (
        mixture.p(),
        p_rgh,
        pimple.dict(),
        false
    ),

    pMin
    (
        "pMin",
        dimPressure,
        mixture
    ),

    rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
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
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        phi*fvc::interpolate(alpha1)
    ),

    K("K", 0.5*magSqr(U)),

    turbulence
    (
        rho,
        U,
        phi,
        rhoPhi,
        alphaPhi1,
        mixture
    ),

    phaseChangePtr
    (
        compressible::twoPhaseChangeModel::New(mixture)
    ),

    phaseChange(*phaseChangePtr),

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
                runTime.timeName(),
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

    if (correctPhi)
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::compressibleVoF::~compressibleVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::compressibleVoF::maxDeltaT() const
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


void Foam::solvers::compressibleVoF::prePredictor()
{
    fvModels().correct();
    alphaPredictor();
    turbulence.correctPhasePhi();
}


void Foam::solvers::compressibleVoF::preSolve()
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

    // Store divU from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    if (correctPhi && mesh.topoChanged())
    {
        // Construct and register divU for mapping
        divU = new volScalarField
        (
            "divU0",
            fvc::div(fvc::absolute(phi, U))
        );
    }

    fvModels().preUpdateMesh();

    // Update the mesh for topology change, mesh to mesh mapping
    mesh.update();
}


void Foam::solvers::compressibleVoF::momentumTransportCorrector()
{
    if (pimple.transportCorr())
    {
        turbulence.correct();
    }
}


void Foam::solvers::compressibleVoF::thermophysicalTransportCorrector()
{}


void Foam::solvers::compressibleVoF::postSolve()
{
    divU.clear();
}


// ************************************************************************* //
