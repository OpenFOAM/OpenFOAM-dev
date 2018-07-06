/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    potentialFoam

Description
    Potential flow solver which solves for the velocity potential, to
    calculate the flux-field, from which the velocity field is obtained by
    reconstructing the flux.

    This application is particularly useful to generate starting fields for
    Navier-Stokes codes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "pName",
        "pName",
        "Name of the pressure field"
    );

    argList::addBoolOption
    (
        "initialiseUBCs",
        "Initialise U boundary conditions"
    );

    argList::addBoolOption
    (
        "writePhi",
        "Write the velocity potential field"
    );

    argList::addBoolOption
    (
        "writep",
        "Calculate and write the pressure field"
    );

    argList::addBoolOption
    (
        "withFunctionObjects",
        "execute functionObjects"
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl potentialFlow(mesh, "potentialFlow");

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Calculating potential flow" << endl;

    // Since solver contains no time loop it would never execute
    // function objects so do it ourselves
    runTime.functionObjects().start();

    MRF.makeRelative(phi);
    adjustPhi(phi, U, p);

    // Non-orthogonal velocity potential corrector loop
    while (potentialFlow.correctNonOrthogonal())
    {
        fvScalarMatrix PhiEqn
        (
            fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
         ==
            fvc::div(phi)
        );

        PhiEqn.setReference(PhiRefCell, PhiRefValue);
        PhiEqn.solve();

        if (potentialFlow.finalNonOrthogonalIter())
        {
            phi -= PhiEqn.flux();
        }
    }

    MRF.makeAbsolute(phi);

    Info<< "Continuity error = "
        << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
        << endl;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

    Info<< "Interpolated velocity error = "
        << (sqrt(sum(sqr(fvc::flux(U) - phi)))/sum(mesh.magSf())).value()
        << endl;

    // Write U and phi
    U.write();
    phi.write();

    // Optionally write Phi
    if (args.optionFound("writePhi"))
    {
        Phi.write();
    }

    // Calculate the pressure field
    if (args.optionFound("writep"))
    {
        Info<< nl << "Calculating approximate pressure field" << endl;

        label pRefCell = 0;
        scalar pRefValue = 0.0;
        setRefCell
        (
            p,
            potentialFlow.dict(),
            pRefCell,
            pRefValue
        );

        // Calculate the flow-direction filter tensor
        volScalarField magSqrU(magSqr(U));
        volSymmTensorField F(sqr(U)/(magSqrU + small*average(magSqrU)));

        // Calculate the divergence of the flow-direction filtered div(U*U)
        // Filtering with the flow-direction generates a more reasonable
        // pressure distribution in regions of high velocity gradient in the
        // direction of the flow
        volScalarField divDivUU
        (
            fvc::div
            (
                F & fvc::div(phi, U),
                "div(div(phi,U))"
            )
        );

        // Solve a Poisson equation for the approximate pressure
        while (potentialFlow.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(p) + divDivUU
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();
        }

        p.write();
    }

    runTime.functionObjects().end();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
