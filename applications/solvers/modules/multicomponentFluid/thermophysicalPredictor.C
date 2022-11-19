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

#include "multicomponentFluid.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multicomponentFluid::thermophysicalPredictor()
{
    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.schemes().div("div(phi,Yi_h)")
        )
    );

    reaction->correct();

    forAll(Y, i)
    {
        if (composition.solve(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              + thermophysicalTransport->divj(Yi)
             ==
                reaction->R(Yi)
              + fvModels().source(rho, Yi)
            );

            YiEqn.relax();

            fvConstraints().constrain(YiEqn);

            YiEqn.solve("Yi");

            fvConstraints().constrain(Yi);
        }
    }

    composition.normalise();


    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + pressureWork
        (
            he.name() == "e"
          ? mvConvection->fvcDiv(phi, p/rho)
          : -dpdt
        )
      + thermophysicalTransport->divq(he)
     ==
        reaction->Qdot()
      + (
            buoyancy.valid()
          ? fvModels().source(rho, he) + rho*(U & buoyancy->g)
          : fvModels().source(rho, he)
        )
    );

    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();

    fvConstraints().constrain(he);

    thermo.correct();
}


// ************************************************************************* //
