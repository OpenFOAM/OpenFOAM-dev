/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "MULES.H"
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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseSystem::solve
(
    const PtrList<volScalarField>& rAUs,
    const PtrList<surfaceScalarField>& rAUfs
)
{
    const dictionary& alphaControls = mesh_.solverDict("alpha");

    const label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));
    const label nAlphaCorr(alphaControls.lookup<label>("nAlphaCorr"));

    const bool LTS = fv::localEulerDdt::enabled(mesh_);

    // Temporary switch for testing and comparing the standard split
    // and the new un-split phase flux discretisation
    const bool splitPhaseFlux
    (
        alphaControls.lookupOrDefault<Switch>("splitPhaseFlux", false)
    );

    // Temporary switch for testing and comparing the standard mean flux
    // and the new phase flux reference for the phase flux correction
    const bool meanFluxReference
    (
        alphaControls.lookupOrDefault<Switch>("meanFluxReference", false)
    );

    // Optional reference phase which is not solved for
    // but obtained from the sum of the other phases
    phaseModel* referencePhasePtr = nullptr;

    // The phases which are solved
    // i.e. the moving phases less the optional reference phase
    phaseModelPartialList solvePhases;

    if (referencePhaseName_ != word::null)
    {
        referencePhasePtr = &phases()[referencePhaseName_];

        solvePhases.setSize(movingPhases().size() - 1);
        label solvePhasesi = 0;
        forAll(movingPhases(), movingPhasei)
        {
            if (&movingPhases()[movingPhasei] != referencePhasePtr)
            {
                solvePhases.set(solvePhasesi++, &movingPhases()[movingPhasei]);
            }
        }
    }
    else
    {
        solvePhases = movingPhases();
    }

    forAll(phases(), phasei)
    {
        phases()[phasei].correctBoundaryConditions();
    }

    PtrList<surfaceScalarField> alphaPhiDbyA0s(phases().size());
    if (implicitPhasePressure() && (rAUs.size() || rAUfs.size()))
    {
        const PtrList<surfaceScalarField> DByAfs(this->DByAfs(rAUs, rAUfs));

        forAll(solvePhases, solvePhasei)
        {
            const phaseModel& phase = solvePhases[solvePhasei];
            const volScalarField& alpha = phase;

            alphaPhiDbyA0s.set
            (
                phase.index(),
                DByAfs[phase.index()]
               *fvc::snGrad(alpha, "bounded")*mesh_.magSf()
            );
        }
    }

    // Calculate the void fraction
    volScalarField alphaVoid
    (
        IOobject
        (
            "alphaVoid",
            mesh_.time().timeName(),
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

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
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
                        "Sp",
                        mesh_.time().timeName(),
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
                    "Su",
                    min(alpha.v(), scalar(1))
                   *fvc::div(fvc::absolute(phi_, phase.U()))->v()
                )
            );

            if (dilatation)
            {
                // Construct the dilatation rate source term
                volScalarField::Internal dgdt
                (
                    volScalarField::Internal::New
                    (
                        "dgdt",
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
                        if (phase.divU().valid())
                        {
                            dgdt += alpha2()*phase.divU()()();
                        }

                        if (phase2.divU().valid())
                        {
                            dgdt -= alpha()*phase2.divU()()();
                        }
                    }
                }

                volScalarField::Internal& Sp = Sps[phasei];
                volScalarField::Internal& Su = Sus[phasei];

                forAll(dgdt, celli)
                {
                    if (dgdt[celli] > 0)
                    {
                        Sp[celli] -= dgdt[celli]/max(1 - alpha[celli], 1e-4);
                        Su[celli] += dgdt[celli]/max(1 - alpha[celli], 1e-4);
                    }
                    else if (dgdt[celli] < 0)
                    {
                        Sp[celli] += dgdt[celli]/max(alpha[celli], 1e-4);
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

        List<volScalarField*> alphaPtrs(phases().size());
        forAll(phases(), phasei)
        {
            alphaPtrs[phasei] = &phases()[phasei];
        }

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
            // Create correction fluxes
            PtrList<surfaceScalarField> alphaPhis(phases().size());

            forAll(movingPhases(), movingPhasei)
            {
                const phaseModel& phase = movingPhases()[movingPhasei];
                const volScalarField& alpha = phase;

                alphaPhis.set
                (
                    phase.index(),
                    new surfaceScalarField
                    (
                        IOobject::groupName("alphaPhiCorr", phase.name()),
                        fvc::flux
                        (
                            splitPhaseFlux ? phi_ : phase.phi()(),
                            alpha,
                            "div(phi," + alpha.name() + ')'
                        )
                    )
                );

                surfaceScalarField& alphaPhi = alphaPhis[phase.index()];

                if (splitPhaseFlux)
                {
                    forAll(phases(), phasei)
                    {
                        const phaseModel& phase2 = phases()[phasei];
                        const volScalarField& alpha2 = phase2;

                        if (&phase2 == &phase) continue;

                        surfaceScalarField phir(phase.phi() - phase2.phi());

                        cAlphaTable::const_iterator cAlpha
                        (
                            cAlphas_.find
                            (
                                phasePairKey(phase.name(), phase2.name())
                            )
                        );

                        if (cAlpha != cAlphas_.end())
                        {
                            surfaceScalarField phic
                            (
                                (mag(phi_) + mag(phir))/mesh_.magSf()
                            );

                            phir +=
                                min(cAlpha()*phic, max(phic))
                               *nHatf(alpha, alpha2);
                        }

                        const word phirScheme
                        (
                            "div(phir,"
                          + alpha2.name() + ',' + alpha.name()
                          + ')'
                        );

                        alphaPhi += fvc::flux
                        (
                           -fvc::flux(-phir, alpha2, phirScheme),
                            alpha,
                            phirScheme
                        );
                    }
                }
                else if (!cAlphas_.empty())
                {
                    forAll(phases(), phasei)
                    {
                        const phaseModel& phase2 = phases()[phasei];
                        const volScalarField& alpha2 = phase2;

                        if (&phase2 == &phase) continue;

                        cAlphaTable::const_iterator cAlpha
                        (
                            cAlphas_.find
                            (
                                phasePairKey(phase.name(), phase2.name())
                            )
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

                if (alphaPhiDbyA0s.set(phase.index()))
                {
                    alphaPhi +=
                        fvc::interpolate(max(alpha, scalar(0)))
                       *fvc::interpolate(max(1 - alpha, scalar(0)))
                       *alphaPhiDbyA0s[phase.index()];
                }

                phase.correctInflowOutflow(alphaPhi);

                MULES::limit
                (
                    geometricOneField(),
                    alpha,
                    meanFluxReference
                      ? phiMoving    // Guarantees boundedness but less accurate
                      : phase.phi()(), // Less robust but more accurate
                    alphaPhi,
                    Sps[phase.index()],
                    Sus[phase.index()],
                    min(alphaVoid.primitiveField(), phase.alphaMax())(),
                    zeroField(),
                    false
                );
            }

            // Limit the flux corrections to ensure the phase fractions sum to 1
            {
                // Generate alphas for the moving phases
                UPtrList<const volScalarField> alphasMoving
                (
                    movingPhases().size()
                );

                UPtrList<surfaceScalarField> alphaPhisMoving
                (
                    movingPhases().size()
                );

                forAll(movingPhases(), movingPhasei)
                {
                    const phaseModel& phase = movingPhases()[movingPhasei];

                    alphasMoving.set(movingPhasei, &phase);

                    alphaPhisMoving.set
                    (
                        movingPhasei,
                        &alphaPhis[phase.index()]
                    );
                }

                MULES::limitSum(alphasMoving, alphaPhisMoving, phiMoving);
            }

            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];
                volScalarField& alpha = phase;

                surfaceScalarField& alphaPhi = alphaPhis[phase.index()];
                phase.correctInflowOutflow(alphaPhi);

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

            if (implicitPhasePressure() && (rAUs.size() || rAUfs.size()))
            {
                const PtrList<surfaceScalarField> DByAfs
                (
                    this->DByAfs(rAUs, rAUfs)
                );

                forAll(solvePhases, solvePhasei)
                {
                    phaseModel& phase = solvePhases[solvePhasei];
                    volScalarField& alpha = phase;

                    const surfaceScalarField alphaDbyA
                    (
                        fvc::interpolate(max(alpha, scalar(0)))
                       *fvc::interpolate(max(1 - alpha, scalar(0)))
                       *DByAfs[phase.index()]
                    );

                    fvScalarMatrix alphaEqn
                    (
                        fvm::ddt(alpha) - fvc::ddt(alpha)
                      - fvm::laplacian(alphaDbyA, alpha, "bounded")
                    );

                    alphaEqn.solve();

                    phase.alphaPhiRef() += alphaEqn.flux();
                }
            }

            // Report the phase fractions and the phase fraction sum
            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];

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
                volScalarField sumAlphaMoving
                (
                    IOobject
                    (
                        "sumAlphaMoving",
                        mesh_.time().timeName(),
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
                    << (sumAlphaMoving + 1 - alphaVoid)()
                      .weightedAverage(mesh_.V()).value()
                    << ' ' << min(sumAlphaMoving + 1 - alphaVoid).value()
                    << ' ' << max(sumAlphaMoving + 1 - alphaVoid).value()
                    << endl;

                // Correct the sum of the phase fractions to avoid drift
                forAll(movingPhases(), movingPhasei)
                {
                    movingPhases()[movingPhasei] *= alphaVoid/sumAlphaMoving;
                }
            }
        }

        if (nAlphaSubCycles > 1)
        {
            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];

                phase.alphaPhiRef() /= nAlphaSubCycles;
            }
        }

        forAll(solvePhases, solvePhasei)
        {
            phaseModel& phase = solvePhases[solvePhasei];

            phase.alphaRhoPhiRef() =
                fvc::interpolate(phase.rho())*phase.alphaPhi();

            phase.maxMin(0, 1);
        }

        if (referencePhasePtr)
        {
            phaseModel& referencePhase = *referencePhasePtr;

            referencePhase.alphaPhiRef() = phi_;

            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];
                referencePhase.alphaPhiRef() -= phase.alphaPhi();
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
