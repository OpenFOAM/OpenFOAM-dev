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

#include "compressibleMultiphaseVoF.H"
#include "subCycle.H"
#include "CMULES.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleMultiphaseVoF::alphaSolve()
{
    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(cAlpha*phic, max(phic));

    UPtrList<const volScalarField> alphas(phases.size());
    PtrList<surfaceScalarField> alphaPhis(phases.size());

    forAll(phases, phasei)
    {
        const compressibleVoFphase& alpha = phases[phasei];

        alphas.set(phasei, &alpha);

        alphaPhis.set
        (
            phasei,
            new surfaceScalarField
            (
                "phi" + alpha.name() + "Corr",
                fvc::flux
                (
                    phi,
                    alpha,
                    alphaScheme
                )
            )
        );

        surfaceScalarField& alphaPhi = alphaPhis[phasei];

        forAll(phases, phasej)
        {
            compressibleVoFphase& alpha2 = phases[phasej];

            if (&alpha2 == &alpha) continue;

            surfaceScalarField phir(phic*mixture.nHatf(alpha, alpha2));

            alphaPhi += fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha,
                alpharScheme
            );
        }

        // Limit alphaPhi for each phase
        MULES::limit
        (
            MULEScontrols,
            1.0/mesh.time().deltaT().value(),
            geometricOneField(),
            alpha,
            phi,
            alphaPhi,
            zeroField(),
            zeroField(),
            oneField(),
            zeroField(),
            false
        );
    }

    MULES::limitSum(alphas, alphaPhis, phi);

    rhoPhi = Zero;

    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    );

    const volScalarField divU(fvc::div(fvc::absolute(phi, U)));

    forAll(phases, phasei)
    {
        compressibleVoFphase& alpha = phases[phasei];

        surfaceScalarField& alphaPhi = alphaPhis[phasei];

        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                mesh.time().name(),
                mesh
            ),
            mesh,
            dimensionedScalar(alpha.vDot().dimensions(), 0)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                mesh.time().name(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            divU.v()*min(alpha.v(), scalar(1))
        );

        {
            const scalarField& vDot = alpha.vDot();

            forAll(vDot, celli)
            {
                if (vDot[celli] < 0.0 && alpha[celli] > 0.0)
                {
                    Sp[celli] += vDot[celli]*alpha[celli];
                    Su[celli] -= vDot[celli]*alpha[celli];
                }
                else if (vDot[celli] > 0.0 && alpha[celli] < 1.0)
                {
                    Sp[celli] -= vDot[celli]*(1.0 - alpha[celli]);
                }
            }
        }


        forAll(phases, phasej)
        {
            const compressibleVoFphase& alpha2 = phases[phasej];

            if (&alpha2 == &alpha) continue;

            const scalarField& vDot2 = alpha2.vDot();

            forAll(vDot2, celli)
            {
                if (vDot2[celli] > 0.0 && alpha2[celli] < 1.0)
                {
                    Sp[celli] -= vDot2[celli]*(1.0 - alpha2[celli]);
                    Su[celli] += vDot2[celli]*alpha[celli];
                }
                else if (vDot2[celli] < 0.0 && alpha2[celli] > 0.0)
                {
                    Sp[celli] += vDot2[celli]*alpha2[celli];
                }
            }
        }

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            alphaPhi,
            Sp,
            Su
        );

        rhoPhi += fvc::interpolate(alpha.thermo().rho())*alphaPhi;

        Info<< alpha.name() << " volume fraction, min, max = "
            << alpha.weightedAverage(mesh.V()).value()
            << ' ' << min(alpha).value()
            << ' ' << max(alpha).value()
            << endl;

        sumAlpha += alpha;
    }

    // Correct the sum of the phase-fractions to avoid 'drift'
    const volScalarField sumCorr(1.0 - sumAlpha);
    forAll(phases, phasei)
    {
        compressibleVoFphase& alpha = phases[phasei];
        alpha += alpha*sumCorr;
    }
}


void Foam::solvers::compressibleMultiphaseVoF::alphaPredictor()
{
    const label nAlphaSubCycles = ceil(nAlphaSubCyclesPtr->value(alphaCoNum));

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.name(),
                mesh
            ),
            mesh,
            dimensionedScalar(rhoPhi.dimensions(), 0)
        );

        const dimensionedScalar totalDeltaT = runTime.deltaT();

        UPtrList<volScalarField> alphas(phases.size());
        forAll(phases, phasei)
        {
            alphas.set(phasei, &phases[phasei]);
        }

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
            alphaSolve();
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        rhoPhi = rhoPhiSum;
    }
    else
    {
        alphaSolve();
    }

    mixture.correct();
}


// ************************************************************************* //
