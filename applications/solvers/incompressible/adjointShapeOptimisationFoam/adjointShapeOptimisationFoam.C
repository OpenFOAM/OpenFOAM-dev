/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    adjointShapeOptimisationFoam

Description
    Steady-state solver for incompressible, turbulent flow of non-Newtonian
    fluids with optimisation of duct shape by applying "blockage" in regions
    causing pressure loss as estimated using an adjoint formulation.

    References:
    \verbatim
        "Implementation of a continuous adjoint for topology optimisation of
         ducted flows"
        C. Othmer,
        E. de Villiers,
        H.G. Weller
        AIAA-2007-3947
        http://pdf.aiaa.org/preview/CDReadyMCFD07_1379/PV2007_3947.pdf
    \endverbatim

    Note that this solver optimises for total pressure loss whereas the
    above paper describes the method for optimising power-loss.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "viscosityModel.H"
#include "incompressibleMomentumTransportModels.H"
#include "simpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"

template<class Type>
void zeroCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = Zero;
    }
}


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
    #include "initAdjointContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        viscosity->lookup("lambda") >> lambda;

        // alpha +=
        //    mesh.relaxationFactor("alpha")
        //   *(lambda*max(Ua & U, zeroSensitivity) - alpha);
        alpha +=
            mesh.fieldRelaxationFactor("alpha")
           *(min(max(alpha + lambda*(Ua & U), zeroAlpha), alphaMax) - alpha);

        zeroCells(alpha, inletCells);
        // zeroCells(alpha, outletCells);

        // Pressure-velocity SIMPLE corrector
        {
            fvModels.correct();

            // Momentum predictor

            tmp<fvVectorMatrix> tUEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevSigma(U)
              + fvm::Sp(alpha, U)
             ==
                fvModels.source(U)
            );
            fvVectorMatrix& UEqn = tUEqn.ref();

            UEqn.relax();

            fvConstraints.constrain(UEqn);

            solve(UEqn == -fvc::grad(p));

            fvConstraints.constrain(U);

            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            tUEqn.clear();
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference
                (
                    pressureReference.refCell(),
                    pressureReference.refValue()
                );

                pEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            fvConstraints.constrain(U);
        }

        // Adjoint Pressure-velocity SIMPLE corrector
        {
            // Adjoint Momentum predictor

            volVectorField adjointTransposeConvection((fvc::grad(Ua) & U));
            // volVectorField adjointTransposeConvection
            //(
            //    fvc::reconstruct
            //    (
            //        mesh.magSf()*fvc::dotInterpolate(fvc::snGrad(Ua), U)
            //    )
            //);

            zeroCells(adjointTransposeConvection, inletCells);

            tmp<fvVectorMatrix> tUaEqn
            (
                fvm::div(-phi, Ua)
              - adjointTransposeConvection
              + turbulence->divDevSigma(Ua)
              + fvm::Sp(alpha, Ua)
             ==
                fvModels.source(Ua)
            );
            fvVectorMatrix& UaEqn = tUaEqn.ref();

            UaEqn.relax();

            fvConstraints.constrain(UaEqn);

            solve(UaEqn == -fvc::grad(pa));

            fvConstraints.constrain(Ua);

            volScalarField rAUa(1.0/UaEqn.A());
            volVectorField HbyAa("HbyAa", Ua);
            HbyAa = rAUa*UaEqn.H();
            tUaEqn.clear();
            surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
            adjustPhi(phiHbyAa, Ua, pa);

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix paEqn
                (
                    fvm::laplacian(rAUa, pa) == fvc::div(phiHbyAa)
                );

                paEqn.setReference(paRefCell, paRefValue);
                paEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phia = phiHbyAa - paEqn.flux();
                }
            }

            #include "adjointContinuityErrs.H"

            // Explicitly relax pressure for adjoint momentum corrector
            pa.relax();

            // Adjoint momentum corrector
            Ua = HbyAa - rAUa*fvc::grad(pa);
            Ua.correctBoundaryConditions();
            fvConstraints.constrain(Ua);
        }

        viscosity->correct();
        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
