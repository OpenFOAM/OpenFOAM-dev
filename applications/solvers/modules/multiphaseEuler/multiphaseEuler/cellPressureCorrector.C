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
    // Face volume fractions
    PtrList<surfaceScalarField> alphafs(phases.size());
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;

        alphafs.set(phasei, fvc::interpolate(alpha).ptr());
        alphafs[phasei].rename("pEqn" + alphafs[phasei].name());
    }

    // Diagonal coefficients
    rAUs.clear();
    rAUs.setSize(phases.size());

    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];
        const volScalarField& alpha = phase;

        rAUs.set
        (
            phase.index(),
            new volScalarField
            (
                IOobject::groupName("rAU", phase.name()),
                1.0
               /(
                   UEqns[phase.index()].A()
                 + byDt
                   (
                       max(phase.residualAlpha() - alpha, scalar(0))
                      *phase.rho()
                   )
                )
            )
        );
    }
    fluid.fillFields("rAU", dimTime/dimDensity, rAUs);

    // Phase diagonal coefficients
    PtrList<surfaceScalarField> alpharAUfs(phases.size());
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;

        alpharAUfs.set
        (
            phasei,
            (
                fvc::interpolate(max(alpha, phase.residualAlpha())*rAUs[phasei])
            ).ptr()
        );
    }

    // Explicit force fluxes
    PtrList<surfaceScalarField> phiFs(fluid.phiFs(rAUs));

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
        fluid.correctBoundaryFlux();

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
                phaseModel& phase = phases[phasei];

                phigFs.set
                (
                    phasei,
                    (
                        alpharAUfs[phasei]
                       *(
                           ghSnGradRho
                         - (fvc::interpolate(phase.rho() - rho))
                          *(buoyancy.g & mesh.Sf())
                         - fluid.surfaceTension(phase)*mesh.magSf()
                        )
                    ).ptr()
                );

                if (phiFs.set(phasei))
                {
                    phigFs[phasei] += phiFs[phasei];
                }
            }
        }

        // Predicted velocities and fluxes for each phase
        PtrList<volVectorField> HbyAs(phases.size());
        PtrList<surfaceScalarField> phiHbyAs(phases.size());
        {
            // Correction force fluxes
            PtrList<surfaceScalarField> ddtCorrByAs(fluid.ddtCorrByAs(rAUs));

            forAll(fluid.movingPhases(), movingPhasei)
            {
                phaseModel& phase = fluid.movingPhases()[movingPhasei];
                const volScalarField& alpha = phase;

                HbyAs.set
                (
                    phase.index(),
                    constrainHbyA
                    (
                        rAUs[phase.index()]
                       *(
                            UEqns[phase.index()].H()
                          + byDt
                            (
                                max(phase.residualAlpha() - alpha, scalar(0))
                               *phase.rho()
                            )
                           *phase.U()().oldTime()
                        ),
                        phase.U(),
                        p_rgh
                    )
                );

                phiHbyAs.set
                (
                    phase.index(),
                    new surfaceScalarField
                    (
                        IOobject::groupName("phiHbyA", phase.name()),
                        fvc::flux(HbyAs[phase.index()])
                      - phigFs[phase.index()]
                      - ddtCorrByAs[phase.index()]
                    )
                );
            }
        }
        fluid.fillFields("HbyA", dimVelocity, HbyAs);
        fluid.fillFields("phiHbyA", dimForce/dimDensity/dimVelocity, phiHbyAs);

        // Add explicit drag forces and fluxes
        PtrList<volVectorField> KdUByAs(fluid.KdUByAs(rAUs));
        PtrList<surfaceScalarField> phiKdPhis(fluid.phiKdPhis(rAUs));

        forAll(phases, phasei)
        {
            if (KdUByAs.set(phasei))
            {
                HbyAs[phasei] -= KdUByAs[phasei];
            }

            if (phiKdPhis.set(phasei))
            {
                phiHbyAs[phasei] -= phiKdPhis[phasei];
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
            dimensionedScalar(dimensionSet(-1, 3, 1, 0, 0), 0)
        );

        forAll(phases, phasei)
        {
            rAUf += alphafs[phasei]*alpharAUfs[phasei];
        }

        rAUf = mag(rAUf);

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
                phaseModel& phase = phases[phasei];
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

                pEqn.solve();
            }

            // Correct fluxes and velocities on last non-orthogonal iteration
            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqnIncomp.flux();

                surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

                forAll(fluid.movingPhases(), movingPhasei)
                {
                    phaseModel& phase = fluid.movingPhases()[movingPhasei];

                    phase.phiRef() =
                        phiHbyAs[phase.index()]
                      + alpharAUfs[phase.index()]*mSfGradp;

                    // Set the phase dilatation rate
                    phase.divU(-pEqnComps[phase.index()] & p_rgh);
                }

                // Optionally relax pressure for velocity correction
                p_rgh.relax();

                mSfGradp = pEqnIncomp.flux()/rAUf;

                forAll(fluid.movingPhases(), movingPhasei)
                {
                    phaseModel& phase = fluid.movingPhases()[movingPhasei];

                    phase.URef() =
                        HbyAs[phase.index()]
                      + fvc::reconstruct
                        (
                            alpharAUfs[phase.index()]*mSfGradp
                          - phigFs[phase.index()]
                        );
                }

                if (partialElimination)
                {
                    fluid.partialElimination(rAUs, KdUByAs, alphafs, phiKdPhis);
                }
                else
                {
                    forAll(fluid.movingPhases(), movingPhasei)
                    {
                        phaseModel& phase = fluid.movingPhases()[movingPhasei];

                        MRF.makeRelative(phase.phiRef());
                        fvc::makeRelative(phase.phiRef(), phase.U());
                    }
                }

                forAll(fluid.movingPhases(), movingPhasei)
                {
                    phaseModel& phase = fluid.movingPhases()[movingPhasei];

                    phase.URef().correctBoundaryConditions();
                    phase.correctUf();
                    fvConstraints().constrain(phase.URef());
                }
            }
        }

        // Update and limit the static pressure
        p = p_rgh + rho*buoyancy.gh;
        fvConstraints().constrain(p);

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
            phaseModel& phase = phases[phasei];
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
