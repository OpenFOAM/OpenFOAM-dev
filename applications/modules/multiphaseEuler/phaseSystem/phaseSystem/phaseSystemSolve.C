/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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

#include "phaseSystem.H"

#include "subCycle.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"
#include "fvcSup.H"

#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"

#include "upwind.H"

#include "CMULES.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseSystem::solve
(
    const alphaControl& alphaControls,
    const PtrList<volScalarField>& rAs
)
{
    const bool boundedPredictor =
       !alphaControls.MULES.globalBounds
     && alphaControls.MULES.extremaCoeff == 0;

    MULES::control MULESBD(alphaControls.MULES);
    MULESBD.globalBounds = true;

    const label nAlphaSubCycles = alphaControls.nAlphaSubCycles;
    const bool LTS = fv::localEulerDdt::enabled(mesh_);

    // Optional reference phase which is not solved for
    // but obtained from the sum of the other phases
    phaseModel* referencePhasePtr = nullptr;

    // The phases which are solved
    // i.e. the moving phases less the optional reference phase
    phaseModelPartialList solvePhases;
    labelList solveMovingPhaseIndices;

    if (referencePhaseName_ != word::null)
    {
        referencePhasePtr = &phases()[referencePhaseName_];

        solvePhases.setSize(movingPhases().size() - 1);
        solveMovingPhaseIndices.setSize(solvePhases.size());

        label solvePhasesi = 0;
        forAll(movingPhases(), movingPhasei)
        {
            if (&movingPhases()[movingPhasei] != referencePhasePtr)
            {
                solveMovingPhaseIndices[solvePhasesi] = movingPhasei;
                solvePhases.set(solvePhasesi++, &movingPhases()[movingPhasei]);
            }
        }
    }
    else
    {
        solvePhases = movingPhases();
        solveMovingPhaseIndices = identityMap(solvePhases.size());
    }

    forAll(phases(), phasei)
    {
        phases()[phasei].correctBoundaryConditions();
    }

    // Create sub-list of alphas and phis for the moving phases
    UPtrList<const volScalarField> alphasMoving(movingPhases().size());
    UPtrList<const surfaceScalarField> phisMoving(movingPhases().size());
    forAll(movingPhases(), movingPhasei)
    {
        alphasMoving.set(movingPhasei, &movingPhases()[movingPhasei]);
        phisMoving.set(movingPhasei, &movingPhases()[movingPhasei].phi()());
    }

    // Calculate the void fraction
    volScalarField alphaVoid
    (
        IOobject
        (
            "alphaVoid",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 1)
    );
    forAll(stationaryPhases(), stationaryPhasei)
    {
        alphaVoid -= stationaryPhases()[stationaryPhasei];
    }

    // Calculate the effective flux of the moving phases
    tmp<surfaceScalarField> tphiMoving(phi_);
    if (stationaryPhases().size())
    {
        tphiMoving = phi_/upwind<scalar>(mesh_, phi_).interpolate(alphaVoid);
    }
    const surfaceScalarField& phiMoving = tphiMoving();

    bool dilatation = false;
    forAll(movingPhases(), movingPhasei)
    {
        if (movingPhases()[movingPhasei].divU().valid())
        {
            dilatation = true;
            break;
        }
    }

    // Optional current phase-pressure diffusion fluxes
    PtrList<surfaceScalarField> alphaPhiDByA(movingPhases().size());

    if (implicitPhasePressure() && rAs.size())
    {
        // Cache the phase-pressure diffusion coefficient
        const surfaceScalarField alphaDByAf(this->alphaDByAf(rAs));

        forAll(movingPhases(), movingPhasei)
        {
            const phaseModel& phase = movingPhases()[movingPhasei];
            const volScalarField& alpha = phase;

            alphaPhiDByA.set
            (
                movingPhasei,
                alphaDByAf*fvc::snGrad(alpha, "bounded")*mesh_.magSf()
            );
        }
    }

    for (int acorr=0; acorr<alphaControls.nAlphaCorr; acorr++)
    {
        PtrList<volScalarField::Internal> Sps(phases().size());
        PtrList<volScalarField::Internal> Sus(phases().size());

        forAll(movingPhases(), movingPhasei)
        {
            const phaseModel& phase = movingPhases()[movingPhasei];
            const volScalarField& alpha = phase;
            const label phasei = phase.index();

            Sps.set
            (
                phasei,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        IOobject::groupName("Sp", phase.name()),
                        mesh_.time().name(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar(dimless/dimTime, 0)
                )
            );

            Sus.set
            (
                phasei,
                new volScalarField::Internal
                (
                    IOobject::groupName("Su", phase.name()),
                    min(alpha.v(), scalar(1))
                   *fvc::div(fvc::absolute(phi_, phase.U()))->v()
                )
            );

            if (dilatation)
            {
                // Construct the dilatation rate source term
                volScalarField::Internal vDot
                (
                    volScalarField::Internal::New
                    (
                        "vDot",
                        mesh_,
                        dimensionedScalar(dimless/dimTime, 0)
                    )
                );

                forAll(phases(), phasej)
                {
                    const phaseModel& phase2 = phases()[phasej];
                    const volScalarField& alpha2 = phase2;

                    if (&phase2 != &phase)
                    {
                        if (!phase.stationary() && phase.divU().valid())
                        {
                            vDot += alpha2()*phase.divU()()();
                        }

                        if (!phase2.stationary() && phase2.divU().valid())
                        {
                            vDot -= alpha()*phase2.divU()()();
                        }
                    }
                }

                volScalarField::Internal& Sp = Sps[phasei];
                volScalarField::Internal& Su = Sus[phasei];

                forAll(vDot, celli)
                {
                    if (vDot[celli] > 0)
                    {
                        Sp[celli] -=
                            vDot[celli]
                           /max
                            (
                                1 - alpha[celli],
                                alphaControls.vDotResidualAlpha
                            );
                        Su[celli] +=
                            vDot[celli]
                           /max
                            (
                                1 - alpha[celli],
                                alphaControls.vDotResidualAlpha
                            );
                    }
                    else if (vDot[celli] < 0)
                    {
                        Sp[celli] +=
                            vDot[celli]
                           /max
                            (
                                alpha[celli],
                                alphaControls.vDotResidualAlpha
                            );
                    }
                }
            }
        }

        tmp<volScalarField> trSubDeltaT;

        if (LTS && nAlphaSubCycles > 1)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh_, nAlphaSubCycles);
        }

        UPtrList<volScalarField> alphas(phases().size());
        forAll(phases(), phasei)
        {
            alphas.set(phasei, &phases()[phasei]);
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
            // Create phase volume-fraction fluxes
            PtrList<surfaceScalarField> alphaPhis(movingPhases().size());

            // Create bounded phase volume-fraction fluxes
            PtrList<surfaceScalarField> alphaPhiBDs(movingPhases().size());

            forAll(movingPhases(), movingPhasei)
            {
                const phaseModel& phase = movingPhases()[movingPhasei];
                const volScalarField& alpha = phase;

                alphaPhis.set
                (
                    movingPhasei,
                    new surfaceScalarField
                    (
                        IOobject::groupName("alphaPhiPred", phase.name()),
                        fvc::flux
                        (
                            phase.phi()(),
                            alpha,
                            "div(phi," + alpha.name() + ')'
                        )
                    )
                );

                surfaceScalarField& alphaPhi = alphaPhis[movingPhasei];

                if (!cAlphas_.empty())
                {
                    forAll(phases(), phasei)
                    {
                        const phaseModel& phase2 = phases()[phasei];
                        const volScalarField& alpha2 = phase2;

                        if (&phase2 == &phase) continue;

                        cAlphaTable::const_iterator cAlpha
                        (
                            cAlphas_.find(phaseInterface(phase, phase2))
                        );

                        if (cAlpha != cAlphas_.end())
                        {
                            const surfaceScalarField phir
                            (
                                phase.phi() - phase2.phi()
                            );

                            const surfaceScalarField phic
                            (
                                (mag(phi_) + mag(phir))/mesh_.magSf()
                            );

                            const surfaceScalarField phirc
                            (
                                min(cAlpha()*phic, max(phic))
                               *nHatf(alpha, alpha2)
                            );

                            const word phirScheme
                            (
                                "div(phir,"
                              + alpha2.name() + ',' + alpha.name()
                              + ')'
                            );

                            alphaPhi += fvc::flux
                            (
                                -fvc::flux(-phirc, alpha2, phirScheme),
                                alpha,
                                phirScheme
                            );
                        }
                    }
                }

                // Add the optional phase-pressure diffusion flux
                if (alphaPhiDByA.set(movingPhasei))
                {
                    alphaPhi += alphaPhiDByA[movingPhasei];
                }

                phase.correctInflowOutflow(alphaPhi);

                if (boundedPredictor)
                {
                    alphaPhiBDs.set
                    (
                        movingPhasei,
                        new surfaceScalarField
                        (
                            IOobject::groupName("alphaPhiBD", phase.name()),
                            upwind<scalar>(mesh_, phase.phi()()).flux(alpha)
                        )
                    );

                    const surfaceScalarField::Boundary& alphaPhiBf =
                        alphaPhi.boundaryField();
                    surfaceScalarField::Boundary& alphaPhiBDBf =
                        alphaPhiBDs[movingPhasei].boundaryFieldRef();

                    forAll(alphaPhiBDBf, patchi)
                    {
                        fvsPatchScalarField& alphaPhiBDPf =
                            alphaPhiBDBf[patchi];

                        if (!alphaPhiBDPf.coupled())
                        {
                            alphaPhiBDPf = alphaPhiBf[patchi];
                        }
                    }
                }
            }

            forAll(movingPhases(), movingPhasei)
            {
                const phaseModel& phase = movingPhases()[movingPhasei];
                const volScalarField& alpha = phase;

                MULES::limit
                (
                    boundedPredictor
                      ? MULESBD
                      : alphaControls.MULES,
                    geometricOneField(),
                    alpha,
                    phiMoving,
                    boundedPredictor
                      ? alphaPhiBDs[movingPhasei]
                      : alphaPhis[movingPhasei],
                    Sps[phase.index()],
                    Sus[phase.index()],
                    min(alphaVoid.primitiveField(), phase.alphaMax())(),
                    zeroField(),
                    false
                );
            }

            // Limit the flux of the phases
            // to ensure the phase fractions sum to 1
            MULES::limitSum
            (
                alphasMoving,
                phisMoving,
                boundedPredictor ? alphaPhiBDs : alphaPhis,
                phiMoving
            );

            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];
                volScalarField& alpha = phase;

                const surfaceScalarField& alphaPhi =
                    (boundedPredictor ? alphaPhiBDs : alphaPhis)
                    [solveMovingPhaseIndices[solvePhasei]];

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha,
                    alphaPhi,
                    Sps[phase.index()],
                    Sus[phase.index()]
                );

                if (alphaSubCycle.index() == 1)
                {
                    phase.alphaPhiRef() = alphaPhi;
                }
                else
                {
                    phase.alphaPhiRef() += alphaPhi;
                }
            }

            if (boundedPredictor)
            {
                if (referencePhasePtr)
                {
                    volScalarField& referenceAlpha = *referencePhasePtr;
                    referenceAlpha = alphaVoid;

                    forAll(solvePhases, solvePhasei)
                    {
                        referenceAlpha -= solvePhases[solvePhasei];
                    }
                }

                // Limit the flux of the phases
                // to ensure the phase fractions sum to 1
                MULES::limitSum(alphasMoving, phisMoving, alphaPhis, phiMoving);

                forAll(movingPhases(), movingPhasei)
                {
                    alphaPhis[movingPhasei] -= alphaPhiBDs[movingPhasei];
                }

                forAll(movingPhases(), movingPhasei)
                {
                    const phaseModel& phase = movingPhases()[movingPhasei];
                    const volScalarField& alpha = phase;

                    MULES::limitCorr
                    (
                        alphaControls.MULES,
                        geometricOneField(),
                        alpha,
                        alphaPhiBDs[movingPhasei],
                        alphaPhis[movingPhasei],
                        Sps[phase.index()],
                        min(alphaVoid.primitiveField(), phase.alphaMax())(),
                        zeroField()
                    );
                }

                // Limit the flux of the phases
                // to ensure the phase fractions sum to 1
                MULES::limitSum(phisMoving, alphaPhis, phiMoving);

                forAll(solvePhases, solvePhasei)
                {
                    phaseModel& phase = solvePhases[solvePhasei];
                    volScalarField& alpha = phase;

                    const surfaceScalarField& alphaPhi =
                        alphaPhis[solveMovingPhaseIndices[solvePhasei]];

                    MULES::correct
                    (
                        geometricOneField(),
                        alpha,
                        alphaPhi,
                        Sps[phase.index()]
                    );

                    phase.alphaPhiRef() += alphaPhi;
                }
            }

            // Add the optional implicit phase pressure contribution
            if (implicitPhasePressure() && rAs.size())
            {
                // Cache the phase-pressure diffusion coefficient
                const surfaceScalarField alphaDByAf(this->alphaDByAf(rAs));

                forAll(solvePhases, solvePhasei)
                {
                    phaseModel& phase = solvePhases[solvePhasei];
                    volScalarField& alpha = phase;

                    fvScalarMatrix alphaEqn
                    (
                        fvm::ddt(alpha) - fvc::ddt(alpha)
                      - fvm::laplacian(alphaDByAf, alpha, "bounded")
                    );

                    alphaEqn.solve();

                    phase.alphaPhiRef() += alphaEqn.flux();
                }
            }

            // Report the phase fractions and the phase fraction sum
            forAll(solvePhases, solvePhasei)
            {
                const phaseModel& phase = solvePhases[solvePhasei];

                Info<< phase.name() << " fraction, min, max = "
                    << phase.weightedAverage(mesh_.V()).value()
                    << ' ' << min(phase).value()
                    << ' ' << max(phase).value()
                    << endl;
            }

            if (referencePhasePtr)
            {
                volScalarField& referenceAlpha = *referencePhasePtr;
                referenceAlpha = alphaVoid;

                forAll(solvePhases, solvePhasei)
                {
                    referenceAlpha -= solvePhases[solvePhasei];
                }
            }
            else
            {
                volScalarField::Internal sumAlphaMoving
                (
                    IOobject
                    (
                        "sumAlphaMoving",
                        mesh_.time().name(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar(dimless, 0)
                );

                forAll(movingPhases(), movingPhasei)
                {
                    sumAlphaMoving += movingPhases()[movingPhasei];
                }

                Info<< "Phase-sum volume fraction, min, max = "
                    << (sumAlphaMoving + 1 - alphaVoid())()
                      .weightedAverage(mesh_.V()).value()
                    << ' ' << min(sumAlphaMoving + 1 - alphaVoid()).value()
                    << ' ' << max(sumAlphaMoving + 1 - alphaVoid()).value()
                    << endl;

                if (alphaControls.clip)
                {
                    // Scale moving phase-fractions to sum to alphaVoid
                    forAll(movingPhases(), movingPhasei)
                    {
                        movingPhases()[movingPhasei].internalFieldRef() *=
                            alphaVoid()/sumAlphaMoving;
                    }
                }
            }
        }

        if (nAlphaSubCycles > 1)
        {
            forAll(solvePhases, solvePhasei)
            {
                solvePhases[solvePhasei].alphaPhiRef() /= nAlphaSubCycles;
            }
        }

        forAll(solvePhases, solvePhasei)
        {
            phaseModel& phase = solvePhases[solvePhasei];
            volScalarField& alpha = phase;

            phase.alphaRhoPhiRef() =
                fvc::interpolate(phase.rho())*phase.alphaPhi();

            if (alphaControls.clip)
            {
                // Clip the phase-fractions between 0 and alphaMax
                alpha.maxMin(0, phase.alphaMax());

                if (stationaryPhases().size())
                {
                    alpha = min(alpha, alphaVoid);
                }
            }
        }

        if (referencePhasePtr)
        {
            phaseModel& referencePhase = *referencePhasePtr;

            referencePhase.alphaPhiRef() = phi_;

            forAll(solvePhases, solvePhasei)
            {
                referencePhase.alphaPhiRef() -=
                    solvePhases[solvePhasei].alphaPhi();
            }

            referencePhase.alphaRhoPhiRef() =
                fvc::interpolate(referencePhase.rho())
               *referencePhase.alphaPhi();

            volScalarField& referenceAlpha = referencePhase;
            referenceAlpha = alphaVoid;

            forAll(solvePhases, solvePhasei)
            {
                referenceAlpha -= solvePhases[solvePhasei];
            }
        }
    }
}


// ************************************************************************* //
