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
#include "momentumTransferSystem.H"

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

bool Foam::phaseSystem::implicitPhasePressure() const
{
    forAll(phases(), phasei)
    {
        if
        (
            mesh()
           .solution()
           .solverDict(phases()[phasei].volScalarField::name())
           .lookupOrDefault<Switch>("implicitPhasePressure", false)
        )
        {
            return true;
        }
    }

    return false;
}


void Foam::phaseSystem::solve
(
    const alphaControl& alphaControls,
    const PtrList<volScalarField>& rAs,
    const momentumTransferSystem& mts
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
    label referenceMovingPhaseIndex = -1;

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
            else
            {
                referenceMovingPhaseIndex = movingPhasei;
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

    // Cache a list of phase-fractions
    UPtrList<volScalarField> alphas(phases().size());
    forAll(phases(), phasei)
    {
        alphas.set(phasei, &phases()[phasei]);
    }

    // Optional current phase-pressure diffusion fluxes
    PtrList<surfaceScalarField> alphaPhiDByA(movingPhases().size());

    if (implicitPhasePressure() && rAs.size())
    {
        // Cache the phase-pressure diffusion coefficient
        const surfaceScalarField alphaDByAf(mts.alphaDByAf(rAs));

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
        // Bounded implicit predictor phase-fractions
        PtrList<volScalarField> alphaPreds(movingPhases().size());

        // Bounded implicit predictor phase-fluxes
        PtrList<surfaceScalarField> alphaPhiPreds(movingPhases().size());

        if (alphaControls.MULESCorr)
        {
            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];
                volScalarField& alpha = phase;

                // Bounded implicit predictor matrix using the mean-flux only
                // to ensure boundedness and phase-consistency
                fvScalarMatrix alphaEqn
                (
                    (
                        LTS
                      ? fv::localEulerDdtScheme<scalar>(mesh_).fvmDdt(alpha)
                      : fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(alpha)
                    )
                  + fv::gaussConvectionScheme<scalar>
                    (
                        mesh_,
                        phiMoving,
                        upwind<scalar>(mesh_, phiMoving)
                    ).fvmDiv(phiMoving, alpha)
                  ==
                    Sus[phase.index()]
                  + fvm::Sp(Sps[phase.index()], alpha)
                );

                alphaEqn.solve();

                alphaPreds.set
                (
                    solveMovingPhaseIndices[solvePhasei],
                    alpha
                );

                alphaPhiPreds.set
                (
                    solveMovingPhaseIndices[solvePhasei],
                    alphaEqn.flux()
                );

                Info<< phase.name() << " predicted fraction, min, max = "
                    << phase.weightedAverage(mesh_.Vsc()).value()
                    << ' ' << min(phase).value()
                    << ' ' << max(phase).value()
                    << endl;
            }

            // Calculate the optional reference phase fraction and flux
            // from the sum of the other phases
            if (referencePhasePtr)
            {
                volScalarField& referenceAlpha = *referencePhasePtr;
                referenceAlpha = alphaVoid;

                forAll(solvePhases, solvePhasei)
                {
                    referenceAlpha -= solvePhases[solvePhasei];
                }

                alphaPreds.set
                (
                    referenceMovingPhaseIndex,
                    referenceAlpha
                );

                alphaPhiPreds.set(referenceMovingPhaseIndex, phi_);

                forAll(solvePhases, solvePhasei)
                {
                    alphaPhiPreds[referenceMovingPhaseIndex]
                        -= alphaPhiPreds[solveMovingPhaseIndices[solvePhasei]];
                }
            }

            // Cache the initial predicted phase-fluxes in the phase
            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];
                const label movingPhasei = solveMovingPhaseIndices[solvePhasei];

                if (alphaSubCycle.index() == 1)
                {
                    phase.alphaPhiRef() = alphaPhiPreds[movingPhasei];
                }
                else
                {
                    phase.alphaPhiRef() += alphaPhiPreds[movingPhasei];
                }
            }
        }

        // Corrector loop
        for (int acorr=0; acorr<alphaControls.nAlphaCorr; acorr++)
        {
            // High-order transport and optional compression phase-fluxes
            PtrList<surfaceScalarField> alphaPhis(movingPhases().size());

            // Bounded phase-fluxes
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
                        IOobject::groupName("alphaPhi", phase.name()),
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

                // Correct the inflow and outflow boundary conditions
                // of the phase flux
                phase.correctInflowOutflow(alphaPhi);

                // Construct the optional bounded phase-flux
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

                    // For non-coupled boundaries
                    // set the bounded flux to the high-order flux
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

            if (alphaControls.MULESCorr)
            {
                // Set local reference to the bounded or high-order phase-fluxes
                PtrList<surfaceScalarField>& alphaPhiCorrs =
                    boundedPredictor ? alphaPhiBDs : alphaPhis;

                forAll(movingPhases(), movingPhasei)
                {
                    const phaseModel& phase = movingPhases()[movingPhasei];
                    const volScalarField& alpha = phase;

                    // Calculate the phase-flux correction
                    // with respect to the bounded implicit prediction
                    alphaPhiCorrs[movingPhasei] -= alphaPhiPreds[movingPhasei];

                    // Limit the phase-flux correction
                    MULES::limitCorr
                    (
                        boundedPredictor
                          ? MULESBD
                          : alphaControls.MULES,
                        geometricOneField(),
                        alpha,
                        alphaPhiPreds[movingPhasei],
                        alphaPhiCorrs[movingPhasei],
                        Sps[phase.index()],
                        min(alphaVoid.primitiveField(), phase.alphaMax())(),
                        zeroField()
                    );
                }

                // Limit the flux corrections of the phases
                // to ensure the phase fractions sum to 1
                MULES::limitSumCorr(alphasMoving, alphaPhiCorrs, phiMoving);

                // Correct the phase-fractions from the phase-flux corrections
                forAll(solvePhases, solvePhasei)
                {
                    phaseModel& phase = solvePhases[solvePhasei];
                    volScalarField& alpha = phase;
                    const label movingPhasei =
                        solveMovingPhaseIndices[solvePhasei];

                    const surfaceScalarField& alphaPhiCorr =
                        alphaPhiCorrs[movingPhasei];

                    MULES::correct
                    (
                        geometricOneField(),
                        alpha,
                        alphaPhiCorr,
                        Sps[phase.index()]
                    );
                }
            }
            else
            {
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
                }
            }

            // Update the phase-fraction of the reference phase
            if (referencePhasePtr)
            {
                volScalarField& referenceAlpha = *referencePhasePtr;
                referenceAlpha = alphaVoid;

                forAll(solvePhases, solvePhasei)
                {
                    referenceAlpha -= solvePhases[solvePhasei];
                }
            }

            if (boundedPredictor)
            {
                // Limit the high-order phase-fluxes
                // to ensure the phase fractions sum to 1
                MULES::limitSum(alphasMoving, alphaPhis, phiMoving);

                // Calculate the phase-flux corrections
                // with respect to the current prediction
                if (alphaControls.MULESCorr)
                {
                    forAll(movingPhases(), movingPhasei)
                    {
                        alphaPhis[movingPhasei] -=
                        (
                            alphaPhiPreds[movingPhasei]
                          + alphaPhiBDs[movingPhasei]
                        );
                    }
                }
                else
                {
                    forAll(movingPhases(), movingPhasei)
                    {
                        alphaPhis[movingPhasei] -= alphaPhiBDs[movingPhasei];
                    }
                }

                // Limit the phase-flux corrections
                forAll(movingPhases(), movingPhasei)
                {
                    const phaseModel& phase = movingPhases()[movingPhasei];
                    const volScalarField& alpha = phase;

                    MULES::limitCorr
                    (
                        alphaControls.MULES,
                        geometricOneField(),
                        alpha,
                        alphaControls.MULESCorr
                          ? alphaPhiPreds[movingPhasei]
                          : alphaPhiBDs[movingPhasei],
                        alphaPhis[movingPhasei],
                        Sps[phase.index()],
                        min(alphaVoid.primitiveField(), phase.alphaMax())(),
                        zeroField()
                    );
                }

                // Limit the flux corrections of the phases
                // to ensure the phase fractions sum to 1
                MULES::limitSumCorr(alphasMoving, alphaPhis, phiMoving);

                // Correct the phase-fractions from the phase-flux corrections
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
                }
            }

            // Add the optional implicit phase pressure contribution
            if (implicitPhasePressure() && rAs.size())
            {
                // Cache the phase-pressure diffusion coefficient
                const surfaceScalarField alphaDByAf(mts.alphaDByAf(rAs));

                // Add the implicit phase-pressure diffusion contribution
                // to the phase-fractions and phase-fluxes
                forAll(solvePhases, solvePhasei)
                {
                    phaseModel& phase = solvePhases[solvePhasei];
                    volScalarField& alpha = phase;

                    fvScalarMatrix alphaEqn
                    (
                        fvm::ddt(alpha) - fvc::ddt(alpha)
                      - fvm::laplacian(alphaDByAf, alpha, "bounded")
                    );

                    alphaEqn.solve
                    (
                        mesh_.solution().solverDict("alpha")
                       .optionalSubDict("phasePressure")
                    );

                    alphaPhis[solveMovingPhaseIndices[solvePhasei]] +=
                        alphaEqn.flux();
                }
            }

            // Report the phase fractions and the phase fraction sum
            forAll(solvePhases, solvePhasei)
            {
                const phaseModel& phase = solvePhases[solvePhasei];

                Info<< phase.name() << " fraction, min, max = "
                    << phase.weightedAverage(mesh_.Vsc()).value()
                    << ' ' << min(phase).value()
                    << ' ' << max(phase).value()
                    << endl;
            }

            // Calculate the reference phase-fraction from the other phases
            // or re-scale the all the phase fractions to sum to 1 exactly
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
                      .weightedAverage(mesh_.Vsc()).value()
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

            if (alphaControls.MULESCorr)
            {
                // For the semi-implicit algorithm relax the phase-fractions
                // and the phase-flux accumulation
                if (boundedPredictor)
                {
                    forAll(movingPhases(), movingPhasei)
                    {
                        phaseModel& phase = movingPhases()[movingPhasei];
                        volScalarField& alpha = phase;

                        alpha +=
                            (1 - alphaControls.MULESCorrRelax)
                           *(alphaPreds[movingPhasei] - alpha);

                        alphaPreds[movingPhasei] = alpha;

                        const surfaceScalarField alphaPhiInc
                        (
                            alphaControls.MULESCorrRelax
                           *(
                               alphaPhiBDs[movingPhasei]
                             + alphaPhis[movingPhasei]
                            )
                        );
                        phase.alphaPhiRef() += alphaPhiInc;
                        alphaPhiPreds[movingPhasei] += alphaPhiInc;
                    }
                }
                else
                {
                    forAll(movingPhases(), movingPhasei)
                    {
                        phaseModel& phase = movingPhases()[movingPhasei];
                        volScalarField& alpha = phase;

                        alpha +=
                            alphaControls.MULESCorrRelax
                           *(alphaPreds[movingPhasei] - alpha);

                        alphaPreds[movingPhasei] = alpha;

                        const surfaceScalarField alphaPhiInc
                        (
                            alphaControls.MULESCorrRelax*alphaPhis[movingPhasei]
                        );
                        phase.alphaPhiRef() += alphaPhiInc;
                        alphaPhiPreds[movingPhasei] += alphaPhiInc;
                    }
                }
            }
            else if (acorr == alphaControls.nAlphaCorr-1)
            {
                // Accumulate the phase-fluxes
                if (boundedPredictor)
                {
                    forAll(movingPhases(), movingPhasei)
                    {
                        phaseModel& phase = movingPhases()[movingPhasei];

                        if (alphaSubCycle.index() == 1)
                        {
                            phase.alphaPhiRef() = alphaPhiBDs[movingPhasei];
                            phase.alphaPhiRef() += alphaPhis[movingPhasei];
                        }
                        else
                        {
                            phase.alphaPhiRef() += alphaPhiBDs[movingPhasei];
                            phase.alphaPhiRef() += alphaPhis[movingPhasei];
                        }
                    }
                }
                else
                {
                    forAll(movingPhases(), movingPhasei)
                    {
                        phaseModel& phase = movingPhases()[movingPhasei];

                        if (alphaSubCycle.index() == 1)
                        {
                            phase.alphaPhiRef() = alphaPhis[movingPhasei];
                        }
                        else
                        {
                            phase.alphaPhiRef() += alphaPhis[movingPhasei];
                        }
                    }
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

    // Calculate the optional reference phase fraction and flux
    // from the sum of the other phases
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


// ************************************************************************* //
