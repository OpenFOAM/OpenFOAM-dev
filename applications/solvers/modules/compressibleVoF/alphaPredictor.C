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

    const volScalarField& rho1 = mixture.thermo1().rho();
    const volScalarField& rho2 = mixture.thermo2().rho();

    tmp<surfaceScalarField> talphaPhi1(alphaPhi1);

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();

        tmp<volScalarField> trSubDeltaT;

        if (LTS)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
        }

        // Create a temporary alphaPhi1 to accumulate the sub-cycled alphaPhi1
        tmp<surfaceScalarField> talphaPhi1
        (
            surfaceScalarField::New
            (
                "alphaPhi1",
                mesh,
                dimensionedScalar(alphaPhi1.dimensions(), 0)
            )
        );

        // Sub-cycle on both alpha1 and alpha2
        List<volScalarField*> alphaPtrs({&alpha1, &alpha2});

        for
        (
            subCycle<volScalarField, subCycleFields> alphaSubCycle
            (
                alphaPtrs,
                nAlphaSubCycles
            );
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi1;
        }

        alphaPhi1 = talphaPhi1();
    }
    else
    {
        #include "alphaEqn.H"
    }

    mixture.correct();

    alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;
    alphaRhoPhi2 = fvc::interpolate(rho2)*(phi - alphaPhi1);

    rhoPhi = alphaRhoPhi1 + alphaRhoPhi2;

    contErr1 =
    (
        fvc::ddt(alpha1, rho1)()() + fvc::div(alphaRhoPhi1)()()
      - (fvModels().source(alpha1, rho1)&rho1)()
    );

    contErr2 =
    (
        fvc::ddt(alpha2, rho2)()() + fvc::div(alphaRhoPhi2)()()
      - (fvModels().source(alpha2, rho2)&rho2)()
    );
}


// ************************************************************************* //
