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

void Foam::solvers::isothermalFluid::correctBuoyantPressure()
{
    // Local references to the buoyancy parameters
    const volScalarField& gh = buoyancy->gh;
    const surfaceScalarField& ghf = buoyancy->ghf;
    const uniformDimensionedScalarField pRef = buoyancy->pRef;

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

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

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

    const bool adjustMass =
        mesh.schemes().steady() && adjustPhi(phiHbyA, U, p_rgh);

    const surfaceScalarField ghGradRhof(-ghf*fvc::snGrad(rho)*mesh.magSf());

    phiHbyA += rhorAUf*ghGradRhof;

    tmp<fvScalarMatrix> tp_rghEqn;

    if (pimple.transonic())
    {
        const surfaceScalarField phidByPsi
        (
            constrainPhid
            (
                fvc::relative(phiHbyA, rho, U)/fvc::interpolate(rho),
                p_rgh
            )
        );

        const surfaceScalarField phid("phid", fvc::interpolate(psi)*phidByPsi);

        // Subtract the compressible part
        // The resulting flux will be zero for a perfect gas
        phiHbyA -= fvc::interpolate(psi*p_rgh)*phidByPsi;

        if (pimple.consistent())
        {
            const surfaceScalarField gradpf(fvc::snGrad(p_rgh)*mesh.magSf());
            phiHbyA += (rhorAAtUf - rhorAUf)*gradpf;
            HbyA += (rAAtU - rAU)*fvc::reconstruct(gradpf - ghGradRhof);
        }

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p_rgh, rho, U, phiHbyA, rhorAAtUf, MRF);

        fvc::makeRelative(phiHbyA, rho, U);

        fvScalarMatrix p_rghDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
          + fvc::div(phiHbyA) + fvm::div(phid, p_rgh)
         ==
            fvModels().source(psi, p_rgh, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            tp_rghEqn = p_rghDDtEqn - fvm::laplacian(rhorAAtUf, p);
            fvScalarMatrix& p_rghEqn = tp_rghEqn.ref();

            // Relax the pressure equation to ensure diagonal-dominance
            p_rghEqn.relax();

            p_rghEqn.setReference
            (
                pressureReference.refCell(),
                pressureReference.refValue()
            );

            p_rghEqn.solve();
        }
    }
    else
    {
        if (pimple.consistent())
        {
            const surfaceScalarField gradpf(fvc::snGrad(p_rgh)*mesh.magSf());
            phiHbyA += (rhorAAtUf - rhorAUf)*gradpf;
            HbyA += (rAAtU - rAU)*fvc::reconstruct(gradpf - ghGradRhof);
        }

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p_rgh, rho, U, phiHbyA, rhorAAtUf, MRF);

        fvc::makeRelative(phiHbyA, rho, U);

        fvScalarMatrix p_rghDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
          + fvc::div(phiHbyA)
         ==
            fvModels().source(psi, p_rgh, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            tp_rghEqn = p_rghDDtEqn - fvm::laplacian(rhorAAtUf, p_rgh);
            fvScalarMatrix& p_rghEqn = tp_rghEqn.ref();

            p_rghEqn.setReference
            (
                pressureReference.refCell(),
                pressureReference.refValue()
            );

            p_rghEqn.solve();
        }
    }

    const fvScalarMatrix& p_rghEqn = tp_rghEqn();

    phi = phiHbyA + p_rghEqn.flux();

    // Calculate and relax the net pressure-buoyancy force
    netForce.ref().relax
    (
        fvc::reconstruct((ghGradRhof + p_rghEqn.flux()/rhorAAtUf)),
        p_rgh.relaxationFactor()
    );

    // Correct the momentum source with the pressure gradient flux
    // calculated from the relaxed pressure
    U = HbyA + rAAtU*netForce();
    U.correctBoundaryConditions();
    fvConstraints().constrain(U);

    K = 0.5*magSqr(U);

    if (!mesh.schemes().steady())
    {
        p = p_rgh + rho*gh + pRef;

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

    p = p_rgh + rho*gh + pRef;

    if (mesh.schemes().steady())
    {
        if (fvConstraints().constrain(p))
        {
            p_rgh = p - rho*gh - pRef;
            p_rgh.correctBoundaryConditions();
        }
    }

    // For steady compressible closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (adjustMass && !thermo.incompressible())
    {
        p += (initialMass - fvc::domainIntegrate(thermo.rho()))
            /fvc::domainIntegrate(psi);
        p_rgh = p - rho*gh - pRef;
        p_rgh.correctBoundaryConditions();
    }

    // Optionally relax pressure for the thermophysics
    p.relax();

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
