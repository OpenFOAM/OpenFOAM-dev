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

#include "incompressibleDenseParticleFluid.H"
#include "fvcDdt.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvmDiv.H"
#include "fvmDdt.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleDenseParticleFluid::momentumPredictor()
{
    volVectorField& Uc(Uc_);

    tUcEqn =
    (
        fvm::ddt(alphac, Uc) + fvm::div(alphaPhic, Uc)
      - fvm::Sp(fvc::ddt(alphac) + fvc::div(alphaPhic), Uc)
      + momentumTransport->divDevSigma(Uc)
     ==
        fvModels().source(alphac, Uc)
    );
    fvVectorMatrix& UcEqn = tUcEqn.ref();

    UcEqn.relax();

    fvConstraints().constrain(UcEqn);

    if (pimple.momentumPredictor())
    {
        // Face buoyancy force
        const surfaceScalarField Fgf(g & mesh.Sf());

        solve
        (
            UcEqn
         ==
            fvc::reconstruct
            (
                Fgf - fvc::snGrad(p)*mesh.magSf()
              - Dcf()*(phic - phid())
            )
          + Dc()*fvc::reconstruct(phic - phid())
          + Fd() - fvm::Sp(Dc(), Uc)
        );

        fvConstraints().constrain(Uc);
    }
}


// ************************************************************************* //
