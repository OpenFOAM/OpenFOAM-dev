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

#include "incompressibleDenseParticleFluid.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "fvcMeshPhi.H"
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleDenseParticleFluid::correctPressure()
{
    volScalarField& p(p_);
    volVectorField& Uc(Uc_);
    surfaceScalarField& phic(phic_);

    fvVectorMatrix& UcEqn = tUcEqn.ref();

    const volScalarField rAUc(1.0/UcEqn.A());
    const volScalarField r1ADUc(1/(1 + rAUc*Dc()));

    const surfaceScalarField rAUcf(fvc::interpolate(rAUc));
    const surfaceScalarField r1ADUcf(1/(1 + rAUcf*Dcf()));
    const surfaceScalarField rADUcf("Dp", r1ADUcf*rAUcf);

    volVectorField HbyA(constrainHbyA(rAUc*UcEqn.H(), Uc, p));

    surfaceScalarField phiHbyAD
    (
        "phiHbyAD",
        (
            r1ADUcf
           *(
                fvc::flux(HbyA)
              + alphacf*rAUcf*fvc::ddtCorr(Uc, phic, Ucf)
            )
        )
    );

    if (p.needReference())
    {
        fvc::makeRelative(phiHbyAD, Uc);
        adjustPhi(phiHbyAD, Uc, p);
        fvc::makeAbsolute(phiHbyAD, Uc);
    }

    // Face buoyancy force
    const surfaceScalarField Fgf(g & mesh.Sf());

    phiHbyAD += rADUcf*(Fgf + Dcf()*phid());

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, Uc, phiHbyAD, rADUcf);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(alphacf*rADUcf, p)
         ==
            fvc::ddt(alphac)
          + fvc::div(alphacf*phiHbyAD)
        );

        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phic = phiHbyAD - pEqn.flux()/alphacf;

            // Explicitly relax pressure for momentum corrector
            p.relax();

            Uc =
                r1ADUc
               *(
                    HbyA
                  + rAUc
                   *(
                        fvc::reconstruct
                        (
                            Fgf - pEqn.flux()/alphacf/rADUcf
                          - Dcf()*(phic - phid())
                        )
                      + Dc()*fvc::reconstruct(phic - phid())
                      + Fd()
                     )
                );

            Uc.correctBoundaryConditions();
            fvConstraints().constrain(Uc);

            // Correct Ucf if the mesh is moving
            fvc::correctUf(Ucf, Uc, phic);

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phic, Uc);
        }
    }
}


// ************************************************************************* //
