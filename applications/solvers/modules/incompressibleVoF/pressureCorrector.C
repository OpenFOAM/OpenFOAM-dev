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

#include "incompressibleVoF.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "findRefCell.H"
#include "fvcMeshPhi.H"
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleVoF::pressureCorrector()
{
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

        if (p_rgh.needReference())
        {
            fvc::makeRelative(phiHbyA, U);
            adjustPhi(phiHbyA, U, p_rgh);
            fvc::makeAbsolute(phiHbyA, U);
        }

        surfaceScalarField phig
        (
            (
                interface.surfaceTensionForce()
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

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix p_rghEqn
            (
                fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
             == Sp_rgh
            );

            p_rghEqn.setReference
            (
                pressureReference().refCell(),
                getRefCellValue(p_rgh, pressureReference().refCell())
            );

            p_rghEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + p_rghEqn.flux();

                p_rgh.relax();

                U = HbyA
                  + rAU()*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
                U.correctBoundaryConditions();
                fvConstraints().constrain(U);
            }
        }

        continuityErrors();

        // Correct Uf if the mesh is moving
        fvc::correctUf(Uf, U, phi, MRF);

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

        p == p_rgh + rho*buoyancy.gh;

        if (p_rgh.needReference())
        {
            p += dimensionedScalar
            (
                "p",
                p.dimensions(),
                pressureReference().refValue()
              - getRefCellValue(p, pressureReference().refCell())
            );
            p_rgh = p - rho*buoyancy.gh;
        }
    }

    if (!correctPhi)
    {
        rAU.clear();
    }

    tUEqn.clear();
}


// ************************************************************************* //
