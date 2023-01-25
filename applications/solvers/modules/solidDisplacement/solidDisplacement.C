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

#include "solidDisplacement.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvmD2dt2.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(solidDisplacement, 0);
    addToRunTimeSelectionTable(solver, solidDisplacement, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::solidDisplacement::readControls()
{
    solid::readControls();

    nCorr = pimple.dict().lookupOrDefault<int>("nCorrectors", 1);
    convergenceTolerance = pimple.dict().lookupOrDefault<scalar>("D", 0);
    pimple.dict().lookup("compactNormalStress") >> compactNormalStress;
    accFac = pimple.dict().lookupOrDefault<scalar>("accelerationFactor", 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::solidDisplacement::solidDisplacement(fvMesh& mesh)
:
    solid
    (
        mesh,
        autoPtr<solidThermo>(new solidDisplacementThermo(mesh))
    ),

    thermo(refCast<solidDisplacementThermo>(solid::thermo)),

    compactNormalStress(pimple.dict().lookup("compactNormalStress")),

    D
    (
        IOobject
        (
            "D",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    E(thermo.E()),
    nu(thermo.nu()),

    mu(E/(2*(1 + nu))),

    lambda
    (
        thermo.planeStress()
      ? nu*E/((1 + nu)*(1 - nu))
      : nu*E/((1 + nu)*(1 - 2*nu))
    ),

    threeK
    (
        thermo.planeStress()
      ? E/(1 - nu)
      : E/(1 - 2*nu)
    ),

    threeKalpha("threeKalpha", threeK*thermo.alphav()),

    sigmaD
    (
        IOobject
        (
            "sigmaD",
            runTime.name(),
            mesh
        ),
        mu*twoSymm(fvc::grad(D)) + lambda*(I*tr(fvc::grad(D)))
    ),

    divSigmaExp
    (
        IOobject
        (
            "divSigmaExp",
            runTime.name(),
            mesh
        ),
        fvc::div(sigmaD)
      - (
            compactNormalStress
          ? fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
          : fvc::div((2*mu + lambda)*fvc::grad(D), "div(sigmaD)")
        )
    )
{
    mesh.schemes().setFluxRequired(D.name());

    // Read the controls
    readControls();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::solidDisplacement::~solidDisplacement()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::solidDisplacement::prePredictor()
{
    if (thermo.thermalStress())
    {
        solid::prePredictor();
    }
}


void Foam::solvers::solidDisplacement::thermophysicalPredictor()
{
    if (thermo.thermalStress())
    {
        solid::thermophysicalPredictor();
    }
}


void Foam::solvers::solidDisplacement::pressureCorrector()
{
    const volScalarField& rho = thermo.rho();

    int iCorr = 0;
    scalar initialResidual = 0;

    {
        {
            fvVectorMatrix DEqn
            (
                fvm::d2dt2(rho, D)
             ==
                fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
              + divSigmaExp
              + rho*fvModels().d2dt2(D)
            );

            if (thermo.thermalStress())
            {
                DEqn += fvc::grad(threeKalpha*T);
            }

            fvConstraints().constrain(DEqn);

            initialResidual = DEqn.solve().max().initialResidual();

            // For steady-state optionally accelerate the solution
            // by over-relaxing the displacement
            if (mesh.schemes().steady() && accFac > 1)
            {
                D += (accFac - 1)*(D - D.oldTime());
            }

            if (!compactNormalStress)
            {
                divSigmaExp = fvc::div(DEqn.flux());
            }
        }

        const volTensorField gradD(fvc::grad(D));
        sigmaD = mu*twoSymm(gradD) + (lambda*I)*tr(gradD);

        if (compactNormalStress)
        {
            divSigmaExp = fvc::div
            (
                sigmaD - (2*mu + lambda)*gradD,
                "div(sigmaD)"
            );
        }
        else
        {
            divSigmaExp += fvc::div(sigmaD);
        }

    } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);
}


void Foam::solvers::solidDisplacement::postCorrector()
{
    if (thermo.thermalStress())
    {
        solid::postCorrector();
    }
}


void Foam::solvers::solidDisplacement::postSolve()
{
    if (runTime.writeTime())
    {
        volSymmTensorField sigma
        (
            IOobject
            (
                "sigma",
                runTime.name(),
                mesh
            ),
            sigmaD
        );

        if (thermo.thermalStress())
        {
            sigma = sigma - I*(threeKalpha*thermo.T());
        }

        volScalarField sigmaEq
        (
            IOobject
            (
                "sigmaEq",
                runTime.name(),
                mesh
            ),
            sqrt((3.0/2.0)*magSqr(dev(sigma)))
        );

        Info<< "Max sigmaEq = " << max(sigmaEq).value()
            << endl;

        sigma.write();
        sigmaEq.write();
    }
}


// ************************************************************************* //
