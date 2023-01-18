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

#include "shockFluid.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::shockFluid::momentumPredictor()
{
    if (!inviscid)
    {
        muEff.clear();
        tauMC.clear();

        muEff = volScalarField::New("muEff", rho*momentumTransport->nuEff());
        tauMC = new volTensorField
        (
            "tauMC",
            muEff()*dev2(Foam::T(fvc::grad(U)))
        );
    }

    const surfaceVectorField phiUp
    (
        (aphiv_pos()*rhoU_pos() + aphiv_neg()*rhoU_neg())
      + (a_pos()*p_pos() + a_neg()*p_neg())*mesh.Sf()
    );

    solve(fvm::ddt(rhoU) + fvc::div(phiUp));

    U.ref() = rhoU()/rho();
    U.correctBoundaryConditions();

    if (!inviscid)
    {
        solve
        (
            fvm::ddt(rho, U) - fvc::ddt(rho, U)
          - fvm::laplacian(muEff(), U)
          - fvc::div(tauMC())
         ==
            fvModels().source(rho, U)
        );

        rhoU == rho*U;
    }
    else
    {
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
    }

    fvConstraints().constrain(U);
}


// ************************************************************************* //
