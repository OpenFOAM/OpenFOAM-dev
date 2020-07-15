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

    // The phases included in the flux sum limit
    // which is all moving phases if the number of solved phases is > 1
    // otherwise it is just the solved phases
    // as the flux sum limit is not needed in this case
    phaseModelPartialList fluxPhases;
    if (solvePhases.size() == 1)
    {
        fluxPhases = solvePhases;
    }
    else
    {
        fluxPhases = movingPhases();
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
            phaseModel& phase = solvePhases[solvePhasei];
            volScalarField& alpha = phase;

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

    bool dilatation = false;
    forAll(fluxPhases, fluxPhasei)
    {
        if (fluxPhases[fluxPhasei].divU().valid())
        {
            dilatation = true;
            break;
        }
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        PtrList<volScalarField::Internal> Sps(phases().size());
        PtrList<volScalarField::Internal> Sus(phases().size());

        forAll(fluxPhases, fluxPhasei)
        {
            phaseModel& phase = fluxPhases[fluxPhasei];
            volScalarField& alpha = phase;
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
                    min(alpha, scalar(1))
                    *fvc::div(fvc::absolute(phi_, phase.U()))
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
            // Generate face-alphas
            PtrList<surfaceScalarField> alphafs(phases().size());
            if (solvePhases.size() > 1)
            {
                forAll(phases(), phasei)
                {
                    phaseModel& phase = phases()[phasei];
                    alphafs.set
                    (
                        phasei,
                        new surfaceScalarField
                        (
                            IOobject::groupName("alphaf", phase.name()),
                            upwind<scalar>(mesh_, phi_).interpolate(phase)
                        )
                    );
                }
            }

            // Create correction fluxes
            PtrList<surfaceScalarField> alphaPhiCorrs(phases().size());

            if (solvePhases.size() > 1)
            {
                forAll(stationaryPhases(), stationaryPhasei)
                {
                    phaseModel& phase = stationaryPhases()[stationaryPhasei];

                    alphaPhiCorrs.set
                    (
                        phase.index(),
                        new surfaceScalarField
                        (
                            IOobject::groupName("alphaPhiCorr", phase.name()),
                          - upwind<scalar>(mesh_, phi_).flux(phase)
                        )
                    );
                }
            }

            forAll(fluxPhases, fluxPhasei)
            {
                phaseModel& phase = fluxPhases[fluxPhasei];
                volScalarField& alpha = phase;

                alphaPhiCorrs.set
                (
                    phase.index(),
                    new surfaceScalarField
                    (
                        IOobject::groupName("alphaPhiCorr", phase.name()),
                        fvc::flux(phi_, alpha, "div(phi," + alpha.name() + ')')
                    )
                );

                surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phase.index()];

                forAll(phases(), phasei)
                {
                    phaseModel& phase2 = phases()[phasei];
                    volScalarField& alpha2 = phase2;

                    if (&phase2 == &phase) continue;

                    surfaceScalarField phir(phase.phi() - phase2.phi());

                    cAlphaTable::const_iterator cAlpha
                    (
                        cAlphas_.find(phasePairKey(phase.name(), phase2.name()))
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

                    word phirScheme
                    (
                        "div(phir," + alpha2.name() + ',' + alpha.name() + ')'
                    );

                    alphaPhiCorr += fvc::flux
                    (
                       -fvc::flux(-phir, alpha2, phirScheme),
                        alpha,
                        phirScheme
                    );
                }

                if (alphaPhiDbyA0s.set(phase.index()))
                {
                    alphaPhiCorr +=
                        fvc::interpolate(max(alpha, scalar(0)))
                       *fvc::interpolate(max(1 - alpha, scalar(0)))
                       *alphaPhiDbyA0s[phase.index()];
                }

                phase.correctInflowOutflow(alphaPhiCorr);

                MULES::limit
                (
                    geometricOneField(),
                    alpha,
                    phi_,
                    alphaPhiCorr,
                    Sps[phase.index()],
                    Sus[phase.index()],
                    min(alphaVoid.primitiveField(), phase.alphaMax())(),
                    zeroField(),
                    true
                );
            }

            if (solvePhases.size() > 1)
            {
                // Limit the flux sums, fixing those of the stationary phases
                labelHashSet fixedAlphaPhiCorrs;
                forAll(stationaryPhases(), stationaryPhasei)
                {
                    fixedAlphaPhiCorrs.insert
                    (
                        stationaryPhases()[stationaryPhasei].index()
                    );
                }
                MULES::limitSum(alphafs, alphaPhiCorrs, fixedAlphaPhiCorrs);
            }

            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];
                volScalarField& alpha = phase;

                surfaceScalarField& alphaPhi = alphaPhiCorrs[phase.index()];
                alphaPhi += upwind<scalar>(mesh_, phi_).flux(phase);
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

            if (!referencePhasePtr)
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

            referencePhase.correctInflowOutflow(referencePhase.alphaPhiRef());
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
