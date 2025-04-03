/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "twoPhaseSolver.H"
#include "subCycle.H"
#include "CMULES.H"
#include "CrankNicolsonDdtScheme.H"
#include "fvcFlux.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::solvers::twoPhaseSolver::alphaPhi
(
    const surfaceScalarField& phi,
    const volScalarField& alpha
)
{
    return fvc::flux
    (
        phi,
        alpha,
        mesh.schemes().div(divAlphaName)
    );
}


void Foam::solvers::twoPhaseSolver::alphaSolve(const label nAlphaSubCycles)
{
    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    {
        tmp<fv::ddtScheme<scalar>> tddtAlpha
        (
            fv::ddtScheme<scalar>::New
            (
                mesh,
                mesh.schemes().ddt("ddt(alpha)")
            )
        );
        const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();

        if
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
         || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
        )
        {
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
        {
            if (nAlphaSubCycles > 1)
            {
                FatalErrorInFunction
                    << "Sub-cycling is not supported "
                       "with the CrankNicolson ddt scheme"
                    << exit(FatalError);
            }

            if
            (
                alphaRestart
             || mesh.time().timeIndex() > mesh.time().startTimeIndex() + 1
            )
            {
                ocCoeff =
                    refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
                   .ocCoeff();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    // Set the time blending factor, 1 for Euler
    const scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = surfaceScalarField::New
        (
            "phiCN",
            cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime()
        );
    }

    tmp<volScalarField> divU;

    if (divergent())
    {
        divU =
        (
            mesh.moving()
          ? fvc::div(phiCN() + mesh.phi())
          : fvc::div(phiCN())
        );
    }

    tmp<volScalarField::Internal> Su;
    tmp<volScalarField::Internal> Sp;

    alphaSuSp(Su, Sp);

    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
              : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiCN,
                upwind<scalar>(mesh, phiCN)
            ).fvmDiv(phiCN, alpha1)
        );

        if (divU.valid())
        {
            alpha1Eqn -= Su() + fvm::Sp(Sp() + divU(), alpha1);
        }

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        tmp<surfaceScalarField> talphaPhi1UD(alpha1Eqn.flux());
        alphaPhi1 = talphaPhi1UD();

        if (alphaApplyPrevCorr && talphaPhi1Corr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct
            (
                MULEScontrols,
                geometricOneField(),
                alpha1,
                alphaPhi1,
                talphaPhi1Corr0.ref(),
                oneField(),
                zeroField()
            );

            alphaPhi1 += talphaPhi1Corr0();
        }

        // Cache the upwind-flux
        talphaPhi1Corr0 = talphaPhi1UD;

        alpha2 = scalar(1) - alpha1;
        alphaPhi2 = phi - alphaPhi1;

        correctInterface();
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        tmp<volScalarField> talpha1CN(alpha1);

        if (ocCoeff > 0)
        {
            // Preserve the BCs of alpha1 in alpha1CN for interpolation
            talpha1CN = alpha1.clone();
            talpha1CN.ref() ==
                (cnCoeff*alpha1 + (1.0 - cnCoeff)*alpha1.oldTime())();
        }

        // Split operator
        tmp<surfaceScalarField> talphaPhi1Un(alphaPhi(phiCN(), talpha1CN()));

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhi1Corr(talphaPhi1Un() - alphaPhi1);
            volScalarField alpha10("alpha10", alpha1);

            if (divU.valid())
            {
                MULES::correct
                (
                    MULEScontrols,
                    geometricOneField(),
                    alpha1,
                    talphaPhi1Un(),
                    talphaPhi1Corr.ref(),
                    (Sp() + divU())(),
                    oneField(),
                    zeroField()
                );
            }
            else
            {
                MULES::correct
                (
                    MULEScontrols,
                    geometricOneField(),
                    alpha1,
                    talphaPhi1Un(),
                    talphaPhi1Corr.ref(),
                    oneField(),
                    zeroField()
                );
            }

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi1 += talphaPhi1Corr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                alphaPhi1 += 0.5*talphaPhi1Corr();
            }
        }
        else
        {
            alphaPhi1 = talphaPhi1Un;

            if (divU.valid())
            {
                MULES::explicitSolve
                (
                    MULEScontrols,
                    geometricOneField(),
                    alpha1,
                    phiCN,
                    alphaPhi1,
                    Sp(),
                    (Su() + divU()*min(alpha1(), scalar(1)))(),
                    oneField(),
                    zeroField()
                );
            }
            else
            {
                MULES::explicitSolve
                (
                    MULEScontrols,
                    geometricOneField(),
                    alpha1,
                    phiCN,
                    alphaPhi1,
                    oneField(),
                    zeroField()
                );
            }
        }

        alpha2 = scalar(1) - alpha1;
        alphaPhi2 = phi - alphaPhi1;

        // Correct only the mixture interface for the interface compression flux
        correctInterface();
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        talphaPhi1Corr0 = alphaPhi1 - talphaPhi1Corr0;

        // Register alphaPhiCorr0.<phase1> for redistribution
        talphaPhi1Corr0.ref().rename
        (
            IOobject::groupName("alphaPhiCorr0", alpha1.group())
        );
        talphaPhi1Corr0.ref().checkIn();
    }
    else
    {
        talphaPhi1Corr0.clear();
    }

    if
    (
        word(mesh.schemes().ddt("ddt(rho,U)"))
     != fv::EulerDdtScheme<vector>::typeName
     && word(mesh.schemes().ddt("ddt(rho,U)"))
     != fv::localEulerDdtScheme<vector>::typeName
    )
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhi1 =
                (alphaPhi1 - (1.0 - cnCoeff)*alphaPhi1.oldTime())/cnCoeff;
            alphaPhi2 = phi - alphaPhi1;
        }
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}


void Foam::solvers::twoPhaseSolver::alphaPredictor()
{
    const label nAlphaSubCycles = ceil(nAlphaSubCyclesPtr->value(alphaCoNum));

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

        UPtrList<volScalarField> alphas({&alpha1, &alpha2});

        for
        (
            subCycle<volScalarField, subCycleFields> alphaSubCycle
            (
                alphas,
                nAlphaSubCycles
            );
            !(++alphaSubCycle).end();
        )
        {
            alphaSolve(nAlphaSubCycles);
            talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi1;
        }

        alphaPhi1 = talphaPhi1();
        alphaPhi2 = phi - talphaPhi1();
    }
    else
    {
        alphaSolve(nAlphaSubCycles);
    }
}


// ************************************************************************* //
