/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

    volScalarField rho("rho", fluid.rho());

    // Correct p_rgh for consistency with the current density
    p_rgh = p - rho*buoyancy.gh - buoyancy.pRef;

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
    rAs.clear();
    if (fluid.implicitPhasePressure())
    {
        rAs.setSize(phases.size());
    }

    PtrList<volVectorField> HVms(movingPhases.size());
    PtrList<PtrList<volScalarField>> invADVs;
    PtrList<PtrList<surfaceScalarField>> invADVfs;
    {
        PtrList<volScalarField> As(movingPhases.size());
        forAll(movingPhases, movingPhasei)
        {
            const phaseModel& phase = movingPhases[movingPhasei];
            const volScalarField& alpha = phase;

            As.set
            (
                movingPhasei,
                UEqns[phase.index()].A()
              + byDt
                (
                    max(phase.residualAlpha() - alpha, scalar(0))
                   *phase.rho()
                )
            );

            if (fluid.implicitPhasePressure())
            {
                rAs.set
                (
                    phase.index(),
                    new volScalarField
                    (
                        IOobject::groupName("rA", phase.name()),
                        1/As[movingPhasei]
                    )
                );
            }
        }

        fluid.invADVs(As, HVms, invADVs, invADVfs);
    }

    // Explicit force fluxes
    PtrList<surfaceScalarField> alphaByADfs;
    PtrList<surfaceScalarField> FgByADfs;
    {
        PtrList<surfaceScalarField> Ffs(fluid.Fs());

        const surfaceScalarField ghSnGradRho
        (
            "ghSnGradRho",
            buoyancy.ghf*fvc::snGrad(rho)*mesh.magSf()
        );

        PtrList<surfaceScalarField> lalphafs(movingPhases.size());
        PtrList<surfaceScalarField> Fgfs(movingPhases.size());

        forAll(movingPhases, movingPhasei)
        {
            const phaseModel& phase = movingPhases[movingPhasei];
            const volScalarField& alpha = phase;

            lalphafs.set
            (
                movingPhasei,
                fvc::interpolate(max(alpha, phase.residualAlpha()))
            );

            Fgfs.set
            (
                movingPhasei,
                (
                    Ffs[phase.index()]
                  + lalphafs[movingPhasei]
                   *(
                       ghSnGradRho
                     - (fvc::interpolate(phase.rho() - rho))
                      *(buoyancy.g & mesh.Sf())
                     - fluid.surfaceTension(phase)*mesh.magSf()
                    )
                ).ptr()
            );
        }

        alphaByADfs = invADVfs & lalphafs;
        FgByADfs = invADVfs & Fgfs;
    }


    // Mass transfer rates
    PtrList<volScalarField> dmdts(fluid.dmdts());
    PtrList<volScalarField> d2mdtdps(fluid.d2mdtdps());

    // --- Optional momentum predictor
    if (predictMomentum)
    {
        PtrList<volVectorField> HbyADs;
        {
            PtrList<volVectorField> Hs(movingPhases.size());

            forAll(movingPhases, movingPhasei)
            {
                const phaseModel& phase = movingPhases[movingPhasei];
                const volScalarField& alpha = phase;

                Hs.set
                (
                    movingPhasei,
                    UEqns[phase.index()].H()
                  + byDt
                    (
                        max(phase.residualAlpha() - alpha, scalar(0))
                       *phase.rho()
                    )
                   *phase.U()().oldTime()
                );

                if (HVms.set(movingPhasei))
                {
                    Hs[movingPhasei] += HVms[movingPhasei];
                }
            }

            HbyADs = invADVs & Hs;
        }

        forAll(movingPhases, movingPhasei)
        {
            const phaseModel& phase = movingPhases[movingPhasei];
            constrainHbyA(HbyADs[movingPhasei], phase.U(), p_rgh);
        }

        const surfaceScalarField mSfGradp(-mesh.magSf()*fvc::snGrad(p_rgh));

        forAll(movingPhases, movingPhasei)
        {
            phaseModel& phase = movingPhases_[movingPhasei];

            phase.URef() =
                HbyADs[movingPhasei]
              + fvc::reconstruct
                (
                    alphaByADfs[movingPhasei]*mSfGradp
                  - FgByADfs[movingPhasei]
                );

            phase.URef().correctBoundaryConditions();
            fvConstraints().constrain(phase.URef());
        }
    }

    // --- Pressure corrector loop
    while (pimple.correct())
    {
        // Correct fixed-flux BCs to be consistent with the velocity BCs
        fluid_.correctBoundaryFlux();

        PtrList<volVectorField> HbyADs;
        PtrList<surfaceScalarField> phiHbyADs;
        {
            // Predicted velocities and fluxes for each phase
            PtrList<volVectorField> Hs(movingPhases.size());
            PtrList<surfaceScalarField> phiHs(movingPhases.size());

            // Correction force fluxes
            PtrList<surfaceScalarField> ddtCorrs(fluid.ddtCorrs());

            forAll(movingPhases, movingPhasei)
            {
                const phaseModel& phase = movingPhases[movingPhasei];
                const volScalarField& alpha = phase;

                Hs.set
                (
                    movingPhasei,
                    UEqns[phase.index()].H()
                  + byDt
                    (
                        max(phase.residualAlpha() - alpha, scalar(0))
                       *phase.rho()
                    )
                   *phase.U()().oldTime()
                );

                if (HVms.set(movingPhasei))
                {
                    Hs[movingPhasei] += HVms[movingPhasei];
                }

                phiHs.set
                (
                    movingPhasei,
                    fvc::flux(Hs[movingPhasei]) + ddtCorrs[phase.index()]
                );
            }

            HbyADs = invADVs & Hs;
            phiHbyADs = invADVfs & phiHs;
        }

        forAll(movingPhases, movingPhasei)
        {
            const phaseModel& phase = movingPhases[movingPhasei];

            constrainHbyA(HbyADs[movingPhasei], phase.U(), p_rgh);
            constrainPhiHbyA(phiHbyADs[movingPhasei], phase.U(), p_rgh);

            phiHbyADs[movingPhasei] -= FgByADfs[movingPhasei];
        }

        // Total predicted flux
        surfaceScalarField phiHbyA
        (
            IOobject
            (
                "phiHbyA",
                runTime.name(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimVolumetricFlux, 0)
        );

        forAll(movingPhases, movingPhasei)
        {
            const phaseModel& phase = movingPhases[movingPhasei];

            phiHbyA += alphafs[phase.index()]*phiHbyADs[movingPhasei];
        }

        MRF.makeRelative(phiHbyA);
        fvc::makeRelative(phiHbyA, movingPhases[0].U());

        // Pressure "diffusivity"
        surfaceScalarField rAf
        (
            IOobject
            (
                "rAf",
                runTime.name(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 0)
        );

        forAll(movingPhases, movingPhasei)
        {
            const phaseModel& phase = movingPhases[movingPhasei];

            rAf += alphafs[phase.index()]*alphaByADfs[movingPhasei];
        }

        // Update the fixedFluxPressure BCs to ensure flux consistency
        {
            surfaceScalarField::Boundary phib
            (
                surfaceScalarField::Internal::null(),
                phi.boundaryField()
            );
            phib = 0;

            forAll(movingPhases, movingPhasei)
            {
                phaseModel& phase = movingPhases_[movingPhasei];

                phib +=
                    alphafs[phase.index()].boundaryField()
                   *phase.phi()().boundaryField();
            }

            setSnGrad<fixedFluxPressureFvPatchScalarField>
            (
                p_rgh.boundaryFieldRef(),
                (
                    phiHbyA.boundaryField() - phib
                )/(mesh.magSf().boundaryField()*rAf.boundaryField())
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
              - fvm::laplacian(rAf, p_rgh)
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

                surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAf);

                forAll(movingPhases, movingPhasei)
                {
                    phaseModel& phase = movingPhases_[movingPhasei];

                    phase.phiRef() =
                        phiHbyADs[movingPhasei]
                      + alphaByADfs[movingPhasei]*mSfGradp;

                    // Set the phase dilatation rate
                    phase.divU(-pEqnComps[phase.index()] & p_rgh);

                    MRF.makeRelative(phase.phiRef());
                    fvc::makeRelative(phase.phiRef(), phase.U());
                }

                // Optionally relax pressure for velocity correction
                p_rgh.relax();

                mSfGradp = pEqnIncomp.flux()/rAf;

                if (dragCorrection)
                {
                    PtrList<volVectorField> dragCorrs(movingPhases.size());
                    PtrList<surfaceScalarField> dragCorrfs(movingPhases.size());
                    fluid.dragCorrs(dragCorrs, dragCorrfs);

                    PtrList<volVectorField> dragCorrByADs
                    (
                        invADVs & dragCorrs
                    );

                    PtrList<surfaceScalarField> dragCorrByADfs
                    (
                        invADVfs & dragCorrfs
                    );

                    forAll(movingPhases, movingPhasei)
                    {
                        phaseModel& phase = movingPhases_[movingPhasei];

                        phase.URef() =
                            HbyADs[movingPhasei]
                          + fvc::reconstruct
                            (
                                alphaByADfs[movingPhasei]*mSfGradp
                              - FgByADfs[movingPhasei]
                              + dragCorrByADfs[movingPhasei]
                            )
                          - dragCorrByADs[movingPhasei];
                    }
                }
                else
                {
                    forAll(movingPhases, movingPhasei)
                    {
                        phaseModel& phase = movingPhases_[movingPhasei];

                        phase.URef() =
                            HbyADs[movingPhasei]
                          + fvc::reconstruct
                            (
                                alphaByADfs[movingPhasei]*mSfGradp
                              - FgByADfs[movingPhasei]
                            );
                    }
                }

                forAll(movingPhases, movingPhasei)
                {
                    phaseModel& phase = movingPhases_[movingPhasei];

                    phase.URef().correctBoundaryConditions();
                    phase.correctUf();
                    fvConstraints().constrain(phase.URef());
                }
            }
        }

        // Update and limit the static pressure
        p_ = p_rgh + rho*buoyancy.gh + buoyancy.pRef;
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
        p_rgh = p - rho*buoyancy.gh - buoyancy.pRef;

        // Update densities from change in p_rgh
        forAll(phases, phasei)
        {
            phaseModel& phase = phases_[phasei];
            if (!phase.incompressible())
            {
                phase.rho() += phase.fluidThermo().psi()*(p_rgh - p_rgh_0);
            }
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
        p_rgh = p - rho*buoyancy.gh - buoyancy.pRef;
    }

    UEqns.clear();
}


// ************************************************************************* //
