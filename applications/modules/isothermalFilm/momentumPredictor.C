/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
#include "fvmDiv.H"
#include "fvcSnGrad.H"
#include "fvcLaplacian.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::solvers::isothermalFilm::sigma() const
{
    return constrainedField(surfaceTension->sigma());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::isothermalFilm::pbByAlphaRhof() const
{
    return fvc::interpolate
    (
        max(nHat & g, dimensionedScalar(g.dimensions(), 0))*VbyA
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::isothermalFilm::pbByAlphaf() const
{
    return fvc::interpolate(rho)*pbByAlphaRhof();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::isothermalFilm::pbByAlphaGradRhof() const
{
    return pbByAlphaRhof()*fvc::snGrad(rho);
}


Foam::tmp<Foam::volScalarField>
Foam::solvers::isothermalFilm::pc(const volScalarField& sigma) const
{
    return -fvc::laplacian(sigma, delta);
}


Foam::tmp<Foam::volScalarField>
Foam::solvers::isothermalFilm::pe() const
{
    // Update the pressure, mapping from the fluid region as required
    p.correctBoundaryConditions();

    // Add the pressure caused normal momentum sources (e.g., parcels impinging
    // with a normal velocity)
    p.internalFieldRef() +=
        VbyA*(nHat & (fvModels().source(alpha, rho, U) & U));

    return p;
}


void Foam::solvers::isothermalFilm::momentumPredictor()
{
    volVectorField& U = U_;

    // Calculate the surface tension coefficient
    const volScalarField sigma(this->sigma());

    // Get the momentum source and remove any normal components
    fvVectorMatrix alphaRhoUsource(fvModels().source(alpha, rho, U));
    alphaRhoUsource.source() -= nHat*(nHat & alphaRhoUsource.source());

    tUEqn =
    (
        fvm::ddt(alpha, rho, U) + fvm::div(alphaRhoPhi, U)
      - fvm::Sp(contErr(), U)
      + momentumTransport->divDevTau(U)
     ==
        contactForce(sigma)
      + alphaRhoUsource
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    // Thermocapillary force
    if (thermocapillary)
    {
        UEqn -= fvc::grad(sigma)/VbyA;
    }

    UEqn.relax();

    fvConstraints().constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        const surfaceScalarField alphaf(fvc::interpolate(alpha));

        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                constrainedField
                (
                    alphaf
                   *(
                        // Buoyancy force
                        fvc::interpolate(rho)*(g & mesh.Sf())

                      - (
                            // External and capillary pressure
                            fvc::snGrad(pe() + pc(sigma), "snGrad(p)")

                            // Buoyant pressure
                          + pbByAlphaGradRhof()*alphaf
                          + pbByAlphaf()*fvc::snGrad(alpha)
                        )*mesh.magSf()
                    )
                )
            )
        );

        // Remove film-normal component of velocity
        U -= nHat*(nHat & U);

        U.correctBoundaryConditions();
    }
}


// ************************************************************************* //
