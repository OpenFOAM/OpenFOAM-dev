/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "multiphaseEuler.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "findRefCell.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSup.H"
#include "fvcSnGrad.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::solvers::multiphaseEuler::cellPressureCorrector()
{
    volScalarField& p(p_);

    // Face volume fractions
    PtrList<surfaceScalarField> alphafs(phases.size());
    forAll(phases, phasei)
    {
        const phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;

        alphafs.set(phasei, fvc::interpolate(max(alpha, scalar(0))).ptr());
        alphafs[phasei].rename("pEqn" + alphafs[phasei].name());
    }

    // Diagonal coefficients
    rAUs.clear();
    rAUs.setSize(phases.size());

    PtrList<surfaceScalarField> rAUfs(phases.size());

    {
        PtrList<volScalarField> Kds(fluid.Kds());

        forAll(fluid.movingPhases(), movingPhasei)
        {
            const phaseModel& phase = fluid.movingPhases()[movingPhasei];
            const volScalarField& alpha = phase;

            volScalarField AU
            (
                UEqns[phase.index()].A()
              + byDt
                (
                    max(phase.residualAlpha() - alpha, scalar(0))
                   *phase.rho()
                )
            );

            if (Kds.set(phase.index()))
            {
                AU += Kds[phase.index()];
            }

            rAUs.set
            (
                phase.index(),
                new volScalarField
                (
                    IOobject::groupName("rAU", phase.name()),
                    1/AU
                )
            );

            rAUfs.set(phase.index(), 1/fvc::interpolate(AU));
        }
        fluid.fillFields("rAU", dimTime/dimDensity, rAUs);
        fluid.fillFields("rAUf", dimTime/dimDensity, rAUfs);
    }

    // Phase diagonal coefficients
    PtrList<surfaceScalarField> alpharAUfs(phases.size());
    forAll(phases, phasei)
    {
        const phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;

        alpharAUfs.set
        (
            phasei,
            (
                fvc::interpolate(max(alpha, phase.residualAlpha()))
               *rAUfs[phasei]
            ).ptr()
        );
    }

    // Explicit force fluxes
    PtrList<surfaceScalarField> Fs(fluid.Fs());

    // Mass transfer rates
    PtrList<volScalarField> dmdts(fluid.dmdts());
    PtrList<volScalarField> d2mdtdps(fluid.d2mdtdps());

    // --- Pressure corrector loop
    while (pimple.correct())
    {
        volScalarField rho("rho", fluid.rho());

        // Correct p_rgh for consistency with p and the updated densities
        p_rgh = p - rho*buoyancy.gh;

        // Correct fixed-flux BCs to be consistent with the velocity BCs
        fluid_.correctBoundaryFlux();

        // Combined buoyancy and force fluxes
        PtrList<surfaceScalarField> phigFs(phases.size());
        {
            const surfaceScalarField ghSnGradRho
            (
                "ghSnGradRho",
                buoyancy.ghf*fvc::snGrad(rho)*mesh.magSf()
            );

            forAll(phases, phasei)
            {
                const phaseModel& phase = phases[phasei];

                phigFs.set
                (
                    phasei,
                    (
                        (
                           ghSnGradRho
                         - (fvc::interpolate(phase.rho() - rho))
                          *(buoyancy.g & mesh.Sf())
                         - fluid.surfaceTension(phase)*mesh.magSf()
                        )
                    ).ptr()
                );
            }
        }

        // Predicted velocities and fluxes for each phase
        PtrList<volVectorField> HbyAs(phases.size());
        PtrList<surfaceScalarField> phiHbyAs(phases.size());
        {
            // Correction force fluxes
            PtrList<surfaceScalarField> ddtCorrs(fluid.ddtCorrs());

            forAll(fluid.movingPhases(), movingPhasei)
            {
                const phaseModel& phase = fluid.movingPhases()[movingPhasei];
                const volScalarField& alpha = phase;
                const label phasei = phase.index();

                const volVectorField H
                (
                    constrainH
                    (
                        UEqns[phasei].H()
                      + byDt
                        (
                            max(phase.residualAlpha() - alpha, scalar(0))
                           *phase.rho()
                        )
                       *phase.U()().oldTime(),
                        rAUs[phasei],
                        phase.U(),
                        p_rgh
                    )
                );

                HbyAs.set(phasei, rAUs[phasei]*H);

                phiHbyAs.set
                (
                    phasei,
                    new surfaceScalarField
                    (
                        IOobject::groupName("phiHbyA", phase.name()),
                        rAUfs[phasei]*(fvc::flux(H) + ddtCorrs[phasei])
                      - alpharAUfs[phasei]*phigFs[phasei]
                      - rAUfs[phasei]*Fs[phasei]
                    )
                );
            }
        }
        fluid.fillFields("HbyA", dimVelocity, HbyAs);
        fluid.fillFields("phiHbyA", dimForce/dimDensity/dimVelocity, phiHbyAs);

        // Add explicit drag forces and fluxes
        PtrList<volVectorField> KdUs(fluid.KdUs());
        PtrList<surfaceScalarField> KdPhis(fluid.KdPhis());

        forAll(phases, phasei)
        {
            if (KdUs.set(phasei))
            {
                HbyAs[phasei] -= rAUs[phasei]*KdUs[phasei];
            }

            if (KdPhis.set(phasei))
            {
                phiHbyAs[phasei] -= rAUfs[phasei]*KdPhis[phasei];
            }
        }

        // Total predicted flux
        surfaceScalarField phiHbyA
        (
            IOobject
            (
                "phiHbyA",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimFlux, 0)
        );

        forAll(phases, phasei)
        {
            phiHbyA += alphafs[phasei]*phiHbyAs[phasei];
        }

        MRF.makeRelative(phiHbyA);
        fvc::makeRelative(phiHbyA, fluid.movingPhases()[0].U());

        // Construct pressure "diffusivity"
        surfaceScalarField rAUf
        (
            IOobject
            (
                "rAUf",
                runTime.name(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 0)
        );

        forAll(phases, phasei)
        {
            rAUf += alphafs[phasei]*alpharAUfs[phasei];
        }

        // Update the fixedFluxPressure BCs to ensure flux consistency
        {
            surfaceScalarField::Boundary phib
            (
                surfaceScalarField::Internal::null(),
                phi.boundaryField()
            );
            phib = 0;

            forAll(phases, phasei)
            {
                const phaseModel& phase = phases[phasei];
                phib +=
                    alphafs[phasei].boundaryField()
                   *phase.phi()().boundaryField();
            }

            setSnGrad<fixedFluxPressureFvPatchScalarField>
            (
                p_rgh.boundaryFieldRef(),
                (
                    phiHbyA.boundaryField() - phib
                )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
            );
        }

        // Compressible pressure equations
        PtrList<fvScalarMatrix> pEqnComps(compressibilityEqns(dmdts, d2mdtdps));

        // Cache p prior to solve for density update
        volScalarField p_rgh_0(p_rgh);

        // Iterate over the pressure equation to correct for non-orthogonality
        while (pimple.correctNonOrthogonal())
        {
            // Construct the transport part of the pressure equation
            fvScalarMatrix pEqnIncomp
            (
                fvc::div(phiHbyA)
              - fvm::laplacian(rAUf, p_rgh)
            );

            // Solve
            {
                fvScalarMatrix pEqn(pEqnIncomp);

                forAll(phases, phasei)
                {
                    pEqn += pEqnComps[phasei];
                }

                if (fluid.incompressible())
                {
                    pEqn.setReference
                    (
                        pressureReference.refCell(),
                        pressureReference.refValue()
                    );
                }

                fvConstraints().constrain(pEqn);

                pEqn.solve();
            }

            // Correct fluxes and velocities on last non-orthogonal iteration
            if (pimple.finalNonOrthogonalIter())
            {
                phi_ = phiHbyA + pEqnIncomp.flux();

                surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

                forAll(fluid.movingPhases(), movingPhasei)
                {
                    phaseModel& phase = fluid_.movingPhases()[movingPhasei];

                    phase.phiRef() =
                        phiHbyAs[phase.index()]
                      + alpharAUfs[phase.index()]*mSfGradp;

                    // Set the phase dilatation rate
                    phase.divU(-pEqnComps[phase.index()] & p_rgh);
                }

                // Optionally relax pressure for velocity correction
                p_rgh.relax();

                mSfGradp = pEqnIncomp.flux()/rAUf;

                if (!dragCorrection)
                {
                    forAll(fluid.movingPhases(), movingPhasei)
                    {
                        phaseModel& phase = fluid_.movingPhases()[movingPhasei];
                        const label phasei = phase.index();

                        phase.URef() =
                            HbyAs[phasei]
                          + fvc::reconstruct
                            (
                                alpharAUfs[phasei]*(mSfGradp - phigFs[phasei])
                              - rAUfs[phasei]*Fs[phasei]
                            );
                    }
                }
                else
                {
                    PtrList<volVectorField> dragCorrs(phases.size());
                    PtrList<surfaceScalarField> dragCorrfs(phases.size());
                    fluid.dragCorrs(dragCorrs, dragCorrfs);

                    forAll(fluid.movingPhases(), movingPhasei)
                    {
                        phaseModel& phase = fluid_.movingPhases()[movingPhasei];
                        const label phasei = phase.index();

                        phase.URef() =
                            HbyAs[phasei]
                          + fvc::reconstruct
                            (
                                alpharAUfs[phasei]*(mSfGradp - phigFs[phasei])
                              + rAUfs[phasei]*(dragCorrfs[phasei] - Fs[phasei])
                            )
                        - rAUs[phasei]*dragCorrs[phasei];
                    }
                }

                if (partialElimination)
                {
                    fluid_.partialElimination
                    (
                        rAUs,
                        KdUs,
                        alphafs,
                        rAUfs,
                        KdPhis
                    );
                }
                else
                {
                    forAll(fluid.movingPhases(), movingPhasei)
                    {
                        phaseModel& phase = fluid_.movingPhases()[movingPhasei];

                        MRF.makeRelative(phase.phiRef());
                        fvc::makeRelative(phase.phiRef(), phase.U());
                    }
                }

                forAll(fluid.movingPhases(), movingPhasei)
                {
                    phaseModel& phase = fluid_.movingPhases()[movingPhasei];

                    phase.URef().correctBoundaryConditions();
                    phase.correctUf();
                    fvConstraints().constrain(phase.URef());
                }
            }
        }

        // Update and limit the static pressure
        p_ = p_rgh + rho*buoyancy.gh;
        fvConstraints().constrain(p_);

        // Account for static pressure reference
        if (p_rgh.needReference() && fluid.incompressible())
        {
            p += dimensionedScalar
            (
                "p",
                p.dimensions(),
                pressureReference.refValue()
              - getRefCellValue(p, pressureReference.refCell())
            );
        }

        // Limit p_rgh
        p_rgh = p - rho*buoyancy.gh;

        // Update densities from change in p_rgh
        forAll(phases, phasei)
        {
            phaseModel& phase = phases_[phasei];
            phase.thermoRef().rho() += phase.thermo().psi()*(p_rgh - p_rgh_0);
        }

        // Update mass transfer rates for change in p_rgh
        forAll(phases, phasei)
        {
            if (dmdts.set(phasei) && d2mdtdps.set(phasei))
            {
                dmdts[phasei] += d2mdtdps[phasei]*(p_rgh - p_rgh_0);
            }
        }

        // Correct p_rgh for consistency with p and the updated densities
        rho = fluid.rho();
        p_rgh = p - rho*buoyancy.gh;
        p_rgh.correctBoundaryConditions();
    }

    UEqns.clear();

    if (!fluid.implicitPhasePressure())
    {
        rAUs.clear();
    }
}


// ************************************************************************* //
