/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "LduMatrix.H"
#include "diagTensorField.H"

typedef LduMatrix<vector, scalar, scalar> lduVectorMatrix;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNo.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        fvVectorMatrix UEqnp(UEqn == -fvc::grad(p));

        lduVectorMatrix U3Eqnp(mesh);
        U3Eqnp.diag() = UEqnp.diag();
        U3Eqnp.upper() = UEqnp.upper();
        U3Eqnp.lower() = UEqnp.lower();
        U3Eqnp.source() = UEqnp.source();

        UEqnp.addBoundaryDiag(U3Eqnp.diag(), 0);
        UEqnp.addBoundarySource(U3Eqnp.source(), false);

        autoPtr<lduVectorMatrix::solver> U3EqnpSolver =
            lduVectorMatrix::solver::New
            (
                U.name(),
                U3Eqnp,
                dictionary
                (
                    IStringStream
                    (
                        "{"
                        "    solver           PBiCG;"
                        "    preconditioner   DILU;"
                        "    tolerance        (1e-5 1e-5 1);"
                        "    relTol           (0 0 0);"
                        "}"
                    )()
                )
            );

        U3EqnpSolver->solve(U).print(Info);

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            volScalarField rAU = 1.0/UEqn.A();

            U = rAU*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf())
                + fvc::ddtCorr(rAU, U, phi);

            adjustPhi(phi, U, p);

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U -= rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
