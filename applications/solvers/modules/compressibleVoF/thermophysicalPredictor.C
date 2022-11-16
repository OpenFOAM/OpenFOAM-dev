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

#include "compressibleVoF.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleVoF::thermophysicalPredictor()
{
    volScalarField& T = mixture.T();
    const volScalarField& p = mixture.p();
    const volScalarField& alpha2(mixture.alpha2());

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T) + fvm::div(rhoPhi, T) - fvm::Sp(contErr(), T)
      - fvm::laplacian(turbulence.alphaEff(), T)
      + (
             mixture.totalInternalEnergy()
          ?
            fvc::div(fvc::absolute(phi, U), p)()() // - contErr()/rho*p
          + (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))()()
          - (U()&(fvModels().source(rho, U)&U)()) - contErr()*K
          :
            p*fvc::div(fvc::absolute(phi, U))()()
        )
       *(
           alpha1()/mixture.thermo1().Cv()()
         + alpha2()/mixture.thermo2().Cv()()
        )
     ==
        fvModels().source(rho, T)
    );

    TEqn.relax();

    fvConstraints().constrain(TEqn);

    TEqn.solve();

    fvConstraints().constrain(T);

    mixture.correctThermo();
    mixture.correct();
}


// ************************************************************************* //
