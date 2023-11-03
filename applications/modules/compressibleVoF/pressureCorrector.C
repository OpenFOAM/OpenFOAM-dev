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

#include "compressibleVoF.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "fvcMeshPhi.H"
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleVoF::pressureCorrector()
{
    volVectorField& U = U_;
    surfaceScalarField& phi(phi_);

    const volScalarField& rho1 = mixture.rho1();
    const volScalarField& rho2 = mixture.rho2();

    const volScalarField& psi1 = mixture.thermo1().psi();
    const volScalarField& psi2 = mixture.thermo2().psi();

    fvVectorMatrix& UEqn = tUEqn.ref();
    setrAU(UEqn);

    const surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));

    const surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);

    while (pimple.correct())
    {
        const volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::flux(HbyA)
          + fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf)
        );

        MRF.makeRelative(phiHbyA);

        const surfaceScalarField phig
        (
            (
                surfaceTensionForce()
              - buoyancy.ghf*fvc::snGrad(rho)
            )*rAUf*mesh.magSf()
        );

        phiHbyA += phig;

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

        // Cache any sources
        fvScalarMatrix p_rghEqnSource
        (
            fvModels().sourceProxy(alpha1, rho1, p_rgh)/rho1
          + fvModels().sourceProxy(alpha2, rho2, p_rgh)/rho2
        );

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phiHbyA, U);

        tmp<fvScalarMatrix> p_rghEqnComp1;
        tmp<fvScalarMatrix> p_rghEqnComp2;

        if (pimple.transonic())
        {
            const surfaceScalarField rho1f(fvc::interpolate(rho1));
            const surfaceScalarField rho2f(fvc::interpolate(rho2));

            surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi);
            surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi);

            p_rghEqnComp1 =
                (
                    (fvc::ddt(alpha1, rho1) + fvc::div(alphaPhi1*rho1f))/rho1
                  - fvc::ddt(alpha1) - fvc::div(alphaPhi1)
                  + (alpha1/rho1)
                   *correction
                    (
                        psi1*fvm::ddt(p_rgh)
                      + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
                    )
                );

            p_rghEqnComp2 =
                (
                   (fvc::ddt(alpha2, rho2) + fvc::div(alphaPhi2*rho2f))/rho2
                 - fvc::ddt(alpha2) - fvc::div(alphaPhi2)
                 + (alpha2/rho2)
                  *correction
                   (
                       psi2*fvm::ddt(p_rgh)
                     + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
                   )
               );
        }
        else
        {
            const surfaceScalarField rho1f(fvc::interpolate(rho1));
            const surfaceScalarField rho2f(fvc::interpolate(rho2));

            p_rghEqnComp1 =
                (
                    (fvc::ddt(alpha1, rho1) + fvc::div(alphaPhi1*rho1f))/rho1
                  - fvc::ddt(alpha1) - fvc::div(alphaPhi1)
                  + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh))
                );

            p_rghEqnComp2 =
                (
                   (fvc::ddt(alpha2, rho2) + fvc::div(alphaPhi2*rho2f))/rho2
                 - fvc::ddt(alpha2) - fvc::div(alphaPhi2)
                 + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh))
                );
        }

        if (mesh.moving())
        {
            p_rghEqnComp1.ref() += fvc::div(mesh.phi())*alpha1;
            p_rghEqnComp2.ref() += fvc::div(mesh.phi())*alpha2;
        }

        p_rghEqnComp1.ref() *= pos(alpha1);
        p_rghEqnComp2.ref() *= pos(alpha2);

        if (pimple.transonic())
        {
            p_rghEqnComp1.ref().relax();
            p_rghEqnComp2.ref().relax();
        }

        // Cache p_rgh prior to solve for density update
        const volScalarField p_rgh_0(p_rgh);

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix p_rghEqnIncomp
            (
                fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
             == p_rghEqnSource
            );

            {
                fvScalarMatrix p_rghEqn
                (
                    p_rghEqnComp1() + p_rghEqnComp2() + p_rghEqnIncomp
                );

                fvConstraints().constrain(p_rghEqn);

                p_rghEqn.solve();
            }

            if (pimple.finalNonOrthogonalIter())
            {
                vDot =
                (
                    alpha1*(p_rghEqnComp2 & p_rgh)
                  - alpha2*(p_rghEqnComp1 & p_rgh)
                );

                phi = phiHbyA + p_rghEqnIncomp.flux();

                p = p_rgh + rho*buoyancy.gh;
                fvConstraints().constrain(p);
                p_rgh = p - rho*buoyancy.gh;
                p_rgh.correctBoundaryConditions();

                U = HbyA
                  + rAU()*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
                U.correctBoundaryConditions();
                fvConstraints().constrain(U);
            }
        }

        // Update densities from change in p_rgh
        mixture_.thermo1().correctRho(psi1*(p_rgh - p_rgh_0));
        mixture_.thermo2().correctRho(psi2*(p_rgh - p_rgh_0));
        mixture_.correct();

        // Correct p_rgh for consistency with p and the updated densities
        p_rgh = p - rho*buoyancy.gh;
        p_rgh.correctBoundaryConditions();
    }

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, fvc::absolute(phi, U), MRF);

    K = 0.5*magSqr(U);

    clearrAU();
    tUEqn.clear();
}


// ************************************************************************* //
