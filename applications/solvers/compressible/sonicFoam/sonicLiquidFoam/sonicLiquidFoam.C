/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    sonicLiquidFoam

Description
    Transient solver for trans-sonic/supersonic, laminar flow of a
    compressible liquid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "compressibleCourantNo.H"

        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U)
              + fvm::div(phi, U)
              - fvm::laplacian(mu, U)
            );

            solve(UEqn == -fvc::grad(p));

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                volScalarField rAU("rAU", 1.0/UEqn.A());
                surfaceScalarField rhorAUf
                (
                    "rhorAUf",
                    fvc::interpolate(rho*rAU)
                );

                U = rAU*UEqn.H();

                surfaceScalarField phid
                (
                    "phid",
                    psi
                   *(
                       fvc::flux(U)
                     + rhorAUf*fvc::ddtCorr(rho, U, phi)/fvc::interpolate(rho)
                    )
                );

                phi = (rhoO/psi)*phid;

                fvScalarMatrix pEqn
                (
                    fvm::ddt(psi, p)
                  + fvc::div(phi)
                  + fvm::div(phid, p)
                  - fvm::laplacian(rhorAUf, p)
                );

                pEqn.solve();

                phi += pEqn.flux();

                solve(fvm::ddt(rho) + fvc::div(phi));
                #include "compressibleContinuityErrs.H"

                U -= rAU*fvc::grad(p);
                U.correctBoundaryConditions();
            }
        }

        rho = rhoO + psi*p;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
