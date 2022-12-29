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
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "fvcMeshPhi.H"
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcVolumeIntegrate.H"
#include "fvmDiv.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::isothermalFluid::correctPressure()
{
    const volScalarField& psi = thermo.psi();
    rho = thermo.rho();
    rho.relax();

    fvVectorMatrix& UEqn = tUEqn.ref();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi*p);

    const volScalarField rAU("rAU", 1.0/UEqn.A());
    const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    tmp<volScalarField> rAtU
    (
        pimple.consistent()
      ? volScalarField::New("rAtU", 1.0/(1.0/rAU - UEqn.H1()))
      : tmp<volScalarField>(nullptr)
    );

    tmp<surfaceScalarField> rhorAtUf
    (
        pimple.consistent()
      ? surfaceScalarField::New("rhoRAtUf", fvc::interpolate(rho*rAtU()))
      : tmp<surfaceScalarField>(nullptr)
    );

    const volScalarField& rAAtU = pimple.consistent() ? rAtU() : rAU;
    const surfaceScalarField& rhorAAtUf =
        pimple.consistent() ? rhorAtUf() : rhorAUf;

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    if (pimple.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf))
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    bool adjustMass = false;

    if (pimple.transonic())
    {
        const surfaceScalarField phidByPsi
        (
            constrainPhid
            (
                fvc::relative(phiHbyA, rho, U)/fvc::interpolate(rho),
                p
            )
        );

        const surfaceScalarField phid("phid", fvc::interpolate(psi)*phidByPsi);

        // Subtract the compressible part
        // The resulting flux will be zero for a perfect gas
        phiHbyA -= fvc::interpolate(psi*p)*phidByPsi;

        if (pimple.consistent())
        {
            phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(p)*mesh.magSf();
            HbyA += (rAAtU - rAU)*fvc::grad(p);
        }

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, rho, U, phiHbyA, rhorAAtUf, MRF);

        fvc::makeRelative(phiHbyA, rho, U);

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyA) + fvm::div(phid, p)
         ==
            fvModels().source(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, p));

            // Relax the pressure equation to ensure diagonal-dominance
            pEqn.relax();

            pEqn.setReference
            (
                pressureReference.refCell(),
                pressureReference.refValue()
            );

            pEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }
    else
    {
        if (pimple.consistent())
        {
            phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(p)*mesh.magSf();
            HbyA += (rAAtU - rAU)*fvc::grad(p);
        }

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, rho, U, phiHbyA, rhorAAtUf, MRF);

        fvc::makeRelative(phiHbyA, rho, U);

        if (mesh.schemes().steady())
        {
            adjustMass = adjustPhi(phiHbyA, U, p);
        }

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyA)
         ==
            fvModels().source(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, p));

            pEqn.setReference
            (
                pressureReference.refCell(),
                pressureReference.refValue()
            );

            pEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    if (!mesh.schemes().steady())
    {
        const bool constrained = fvConstraints().constrain(p);

        // Thermodynamic density update
        thermo.correctRho(psi*p - psip0);

        if (constrained)
        {
            rho = thermo.rho();
        }

        correctDensity();
    }

    continuityErrors();

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAAtU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvConstraints().constrain(U);
    K = 0.5*magSqr(U);

    if (mesh.schemes().steady())
    {
        fvConstraints().constrain(p);
    }

    // For steady compressible closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (adjustMass && !thermo.incompressible())
    {
        p += (initialMass - fvc::domainIntegrate(thermo.rho()))
            /fvc::domainIntegrate(psi);
        p.correctBoundaryConditions();
    }

    if (mesh.schemes().steady() || pimple.simpleRho() || adjustMass)
    {
        rho = thermo.rho();
    }

    // Correct rhoUf if the mesh is moving
    fvc::correctRhoUf(rhoUf, rho, U, phi, MRF);

    if (mesh.schemes().steady() || pimple.simpleRho())
    {
        rho.relax();
    }

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);
    }
}


// ************************************************************************* //
