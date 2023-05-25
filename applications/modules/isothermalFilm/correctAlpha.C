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

#include "isothermalFilm.H"
#include "surfaceTensionModel.H"
#include "constrainHbyA.H"
#include "fvcFlux.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::isothermalFilm::correctAlpha()
{
    volScalarField& alpha = alpha_;

    fvVectorMatrix& UEqn = tUEqn.ref();

    const surfaceScalarField rhof(fvc::interpolate(rho));

    const surfaceScalarField pbByAlphaf(this->pbByAlphaf());
    const surfaceScalarField pbByAlphaGradRhof
    (
        constrainedField(this->pbByAlphaGradRhof()*mesh.magSf())
    );

    const surfaceScalarField phip
    (
        constrainedField
        (
            fvc::snGrad(pe() + pc(surfaceTension->sigma()), "snGrad(p)")
           *mesh.magSf()
          - rhof*(g & mesh.Sf())
        )
    );

    while (pimple.correct())
    {
        const volScalarField rAU(1/UEqn.A());
        const volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, alpha));

        const surfaceScalarField alphaf
        (
            constrainedField(fvc::interpolate(alpha))
        );

        const surfaceScalarField alpharAUf
        (
            constrainedField(fvc::interpolate(alpha*rAU))
        );

        const surfaceScalarField phig("phig", phip + pbByAlphaGradRhof*alphaf);

        phi_ = constrainedField(fvc::flux(HbyA) - alpharAUf*phig);

        const surfaceScalarField phid("phid", rhof*phi);

        const surfaceScalarField rAUf
        (
            "rAUf",
            alphaf*rhof*alpharAUf*pbByAlphaf
        );

        fvScalarMatrix alphadEqn
        (
            fvm::ddt(rho, alpha)
          + fvm::div(phid, alpha)
         ==
            fvModels().source(rho, alpha)
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix alphaEqn(alphadEqn - fvm::laplacian(rAUf, alpha));

            alphaEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                alphaRhoPhi_ = alphaEqn.flux();
            }
        }

        const surfaceScalarField phiGradAlpha
        (
            constrainedField(pbByAlphaf*fvc::snGrad(alpha)*mesh.magSf())
        );

        phi_ -= alpharAUf*phiGradAlpha;

        U_ = HbyA - rAU*fvc::reconstruct(alphaf*(phig + phiGradAlpha));

        // Remove film-normal component of velocity
        U_ -= nHat*(nHat & U);

        U_.correctBoundaryConditions();

        fvConstraints().constrain(U_);

        fvConstraints().constrain(alpha);

        continuityErrors();
    }

    // Update film thickness
    correctDelta();

    tUEqn.clear();
}


// ************************************************************************* //
