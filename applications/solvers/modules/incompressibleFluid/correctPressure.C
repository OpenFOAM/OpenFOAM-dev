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

#include "incompressibleFluid.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "fvcMeshPhi.H"
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleFluid::correctPressure()
{
    fvVectorMatrix& UEqn = tUEqn.ref();

    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf))
    );

    MRF.makeRelative(phiHbyA);

    if (p.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p);
        fvc::makeAbsolute(phiHbyA, U);
    }

    tmp<volScalarField> rAtU(rAU);

    if (pimple.consistent())
    {
        rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    if (pimple.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU(), MRF);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );

        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    continuityErrors();

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAtU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvConstraints().constrain(U);

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi, MRF);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}


// ************************************************************************* //
