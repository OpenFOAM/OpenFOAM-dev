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
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleVoF::pressureCorrector()
{
    volScalarField& p = mixture.p();

    const volScalarField& alpha2(mixture.alpha2());

    const volScalarField& rho1 = mixture.rho1();
    const volScalarField& rho2 = mixture.rho2();

    const volScalarField& psi1 = mixture.thermo1().psi();
    const volScalarField& psi2 = mixture.thermo2().psi();

    fvVectorMatrix& UEqn = tUEqn.ref();

    if (rAU.valid())
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));

    const surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);

    while (pimple.correct())
    {
        volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::flux(HbyA)
          + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
        );

        MRF.makeRelative(phiHbyA);

        surfaceScalarField phig
        (
            (
                mixture.surfaceTensionForce()
              - buoyancy.ghf*fvc::snGrad(rho)
            )*rAUf*mesh.magSf()
        );

        phiHbyA += phig;

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

        // Cache the phase change pressure source
        fvScalarMatrix Sp_rgh
        (
            fvModels().source
            (
                volScalarField::New
                (
                    "1",
                    mesh,
                    dimensionedScalar(dimless/dimPressure, 1)
                ),
                p_rgh
            )
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

        p_rghEqnComp1.ref() -=
            (fvModels().source(alpha1, mixture.thermo1().rho())&rho1)/rho1;
        p_rghEqnComp2.ref() -=
            (fvModels().source(alpha2, mixture.thermo2().rho())&rho2)/rho2;

        if (pimple.transonic())
        {
            p_rghEqnComp1.ref().relax();
            p_rghEqnComp2.ref().relax();
        }

        // Cache p_rgh prior to solve for density update
        volScalarField p_rgh_0(p_rgh);

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix p_rghEqnIncomp
            (
                fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
             == Sp_rgh
            );

            solve
            (
                p_rghEqnComp1() + p_rghEqnComp2() + p_rghEqnIncomp
            );

            if (pimple.finalNonOrthogonalIter())
            {
                dgdt =
                (
                    alpha1*(p_rghEqnComp2 & p_rgh)
                  - alpha2*(p_rghEqnComp1 & p_rgh)
                );

                phi = phiHbyA + p_rghEqnIncomp.flux();

                p = max(p_rgh + (alpha1*rho1 + alpha2*rho2)*buoyancy.gh, pMin);
                p_rgh = p - (alpha1*rho1 + alpha2*rho2)*buoyancy.gh;
                p_rgh.correctBoundaryConditions();

                U = HbyA
                  + rAU()*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
                U.correctBoundaryConditions();
                fvConstraints().constrain(U);
            }
        }

        // Correct Uf if the mesh is moving
        fvc::correctUf(Uf, U, fvc::absolute(phi, U), MRF);

        // Update densities from change in p_rgh
        mixture.thermo1().correctRho(psi1*(p_rgh - p_rgh_0));
        mixture.thermo2().correctRho(psi2*(p_rgh - p_rgh_0));
        mixture.correct();

        // Correct p_rgh for consistency with p and the updated densities
        p_rgh = p - rho*buoyancy.gh;
        p_rgh.correctBoundaryConditions();
    }

    K = 0.5*magSqr(U);

    tUEqn.clear();
}


// ************************************************************************* //
