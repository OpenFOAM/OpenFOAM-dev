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

#include "compressibleMultiphaseVoF.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "findRefCell.H"
#include "fvcMeshPhi.H"
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcSup.H"
#include "fvcReconstruct.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleMultiphaseVoF::pressureCorrector()
{
    volVectorField& U = U_;
    surfaceScalarField& phi(phi_);

    fvVectorMatrix& UEqn = tUEqn.ref();
    setrAU(UEqn);

    const surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));

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

        PtrList<fvScalarMatrix> p_rghEqnComps(phases.size());

        forAll(phases, phasei)
        {
            const compressibleVoFphase& phase = phases[phasei];
            const rhoThermo& thermo = phase.thermo();
            const volScalarField& rho = phases[phasei].thermo().rho();

            p_rghEqnComps.set
            (
                phasei,
                (
                    fvc::ddt(rho) + thermo.psi()*correction(fvm::ddt(p_rgh))
                  + fvc::div(phi, rho) - fvc::Sp(fvc::div(phi), rho)
                  - (fvModels().source(phase, rho)&rho)
                ).ptr()
            );
        }

        // Cache p_rgh prior to solve for density update
        const volScalarField p_rgh_0(p_rgh);

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix p_rghEqnIncomp
            (
                fvc::div(phiHbyA)
              - fvm::laplacian(rAUf, p_rgh)
            );

            tmp<fvScalarMatrix> p_rghEqnComp;

            forAll(phases, phasei)
            {
                const compressibleVoFphase& phase = phases[phasei];

                tmp<fvScalarMatrix> p_rghEqnCompi
                (
                    (max(phase, scalar(0))/phase.thermo().rho())
                   *p_rghEqnComps[phasei]
                );

                if (phasei == 0)
                {
                    p_rghEqnComp = p_rghEqnCompi;
                }
                else
                {
                    p_rghEqnComp.ref() += p_rghEqnCompi;
                }
            }

            {
                fvScalarMatrix p_rghEqn(p_rghEqnComp + p_rghEqnIncomp);

                fvConstraints().constrain(p_rghEqn);

                p_rghEqn.solve();
            }

            if (pimple.finalNonOrthogonalIter())
            {
                forAll(phases, phasei)
                {
                    compressibleVoFphase& phase = phases[phasei];

                    phase.dgdt() =
                        pos0(phase)
                       *(p_rghEqnComps[phasei] & p_rgh)/phase.thermo().rho();
                }

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
        mixture.correctRho(p_rgh - p_rgh_0);
        mixture.correct();

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
