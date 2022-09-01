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
#include "interfaceCompression.H"
#include "CMULES.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleVoF::alphaPredictor()
{
    #include "alphaControls.H"

    volScalarField& alpha2(mixture.alpha2());

    const volScalarField& rho1 = mixture.thermo1().rho();
    const volScalarField& rho2 = mixture.thermo2().rho();

    tmp<surfaceScalarField> talphaPhi1(alphaPhi1);

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();

        talphaPhi1 = new surfaceScalarField
        (
            IOobject
            (
                "alphaPhi1",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(alphaPhi1.dimensions(), 0)
        );

        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(rhoPhi.dimensions(), 0)
        );

        tmp<volScalarField> trSubDeltaT;

        if (LTS)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
        }

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi1;
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        alphaPhi1 = talphaPhi1();
        rhoPhi = rhoPhiSum;
    }
    else
    {
        #include "alphaEqn.H"
    }

    contErr =
        (
            fvc::ddt(rho)()() + fvc::div(rhoPhi)()()
          - (fvModels().source(alpha1, rho1)&rho1)()
          - (fvModels().source(alpha2, rho2)&rho2)()
        );
}


// ************************************************************************* //
