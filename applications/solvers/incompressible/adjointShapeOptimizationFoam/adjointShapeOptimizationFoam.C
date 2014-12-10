/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    ajointShapeOptimizationFoam

Description
    Steady-state solver for incompressible, turbulent flow of non-Newtonian
    fluids with optimisation of duct shape by applying "blockage" in regions
    causing pressure loss as estimated using an adjoint formulation.

    References:
    \verbatim
        "Implementation of a continuous adjoint for topology optimization of
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
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"

template<class Type>
void zeroCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = pTraits<Type>::zero;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "initAdjointContinuityErrs.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        laminarTransport.lookup("lambda") >> lambda;

        //alpha +=
        //    mesh.relaxationFactor("alpha")
        //   *(lambda*max(Ua & U, zeroSensitivity) - alpha);
        alpha +=
            mesh.fieldRelaxationFactor("alpha")
           *(min(max(alpha + lambda*(Ua & U), zeroAlpha), alphaMax) - alpha);

        zeroCells(alpha, inletCells);
        //zeroCells(alpha, outletCells);

        // Pressure-velocity SIMPLE corrector
        {
            // Momentum predictor

            tmp<fvVectorMatrix> UEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevReff(U)
              + fvm::Sp(alpha, U)
            );

            UEqn().relax();

            solve(UEqn() == -fvc::grad(p));

            volScalarField rAU(1.0/UEqn().A());
            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn().H();
            UEqn.clear();
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::interpolate(HbyA) & mesh.Sf()
            );
            adjustPhi(phiHbyA, U, p);

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);
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
        }

        // Adjoint Pressure-velocity SIMPLE corrector
        {
            // Adjoint Momentum predictor

            volVectorField adjointTransposeConvection((fvc::grad(Ua) & U));
            //volVectorField adjointTransposeConvection
            //(
            //    fvc::reconstruct
            //    (
            //        mesh.magSf()*(fvc::snGrad(Ua) & fvc::interpolate(U))
            //    )
            //);

            zeroCells(adjointTransposeConvection, inletCells);

            tmp<fvVectorMatrix> UaEqn
            (
                fvm::div(-phi, Ua)
              - adjointTransposeConvection
              + turbulence->divDevReff(Ua)
              + fvm::Sp(alpha, Ua)
            );

            UaEqn().relax();

            solve(UaEqn() == -fvc::grad(pa));

            volScalarField rAUa(1.0/UaEqn().A());
            volVectorField HbyAa("HbyAa", Ua);
            HbyAa = rAUa*UaEqn().H();
            UaEqn.clear();
            surfaceScalarField phiHbyAa
            (
                "phiHbyAa",
                fvc::interpolate(HbyAa) & mesh.Sf()
            );
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
        }

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
