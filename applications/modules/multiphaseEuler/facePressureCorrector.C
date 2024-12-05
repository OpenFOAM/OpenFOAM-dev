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

void Foam::solvers::multiphaseEuler::facePressureCorrector()
{
    volScalarField& p(p_);

    // Face volume fractions
    PtrList<surfaceScalarField> alphafs(phases.size());
    forAll(phases, phasei)
    {
        const phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;

        alphafs.set(phasei, fvc::interpolate(alpha).ptr());
        alphafs[phasei].rename("pEqn" + alphafs[phasei].name());
    }

    // Diagonal coefficients
    rAs.clear();
    if (fluid.implicitPhasePressure())
    {
        rAs.setSize(phases.size());
    }

    PtrList<surfaceScalarField> HVmfs(movingPhases.size());
    PtrList<PtrList<surfaceScalarField>> invADVfs;
    {
        PtrList<surfaceScalarField> Afs(movingPhases.size());

        forAll(movingPhases, movingPhasei)
        {
            const phaseModel& phase = movingPhases[movingPhasei];
            const volScalarField& alpha = phase;

            const volScalarField A
            (
                byDt
                (
                    max(alpha.oldTime(), phase.residualAlpha())
                   *phase.rho().oldTime()
                )
              + UEqns[phase.index()].A()
            );

            if (fluid.implicitPhasePressure())
            {
                rAs.set
                (
                    phase.index(),
                    new volScalarField
                    (
                        IOobject::groupName("rA", phase.name()),
                        1/A
                    )
                );
            }

            Afs.set
            (
                movingPhasei,
                new surfaceScalarField
                (
                    IOobject::groupName("rAf", phase.name()),
                    fvc::interpolate(A)
                )
            );
        }

        invADVfs = fluid.invADVfs(Afs, HVmfs);
    }

    volScalarField rho("rho", fluid.rho());

    // Phase diagonal coefficients
    PtrList<surfaceScalarField> alphaByADfs;
    PtrList<surfaceScalarField> FgByADfs;
    {
        // Explicit force fluxes
        PtrList<surfaceScalarField> Ffs(fluid.Ffs());

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

            lalphafs.set
            (
                movingPhasei,
                max(alphafs[phase.index()], phase.residualAlpha())
            );

            Fgfs.set
            (
                movingPhasei,
                Ffs[phase.index()]
              + lalphafs[movingPhasei]
               *(
                    ghSnGradRho
                  - (fvc::interpolate(phase.rho() - rho))
                   *(buoyancy.g & mesh.Sf())
                  - fluid.surfaceTension(phase)*mesh.magSf()
                )
            );
        }

        alphaByADfs = invADVfs & lalphafs;
        FgByADfs = invADVfs & Fgfs;
    }


    // Mass transfer rates
    PtrList<volScalarField> dmdts(fluid.dmdts());
    PtrList<volScalarField> d2mdtdps(fluid.d2mdtdps());

    // --- Pressure corrector loop
    while (pimple.correct())
    {
        // Correct fixed-flux BCs to be consistent with the velocity BCs
        fluid_.correctBoundaryFlux();

        // Predicted fluxes for each phase
        PtrList<surfaceScalarField> phiHbyADs;
        {
            PtrList<surfaceScalarField> phiHs(movingPhases.size());

            forAll(movingPhases, movingPhasei)
            {
                const phaseModel& phase = movingPhases[movingPhasei];
                const volScalarField& alpha = phase;

                phiHs.set
                (
                    movingPhasei,
                    (
                       fvc::interpolate
                       (
                           max(alpha.oldTime(), phase.residualAlpha())
                          *phase.rho().oldTime()
                       )
                      *byDt
                       (
                           phase.Uf().valid()
                         ? (mesh.Sf() & phase.Uf()().oldTime())
                         : MRF.absolute(phase.phi()().oldTime())
                       )
                     + fvc::flux(UEqns[phase.index()].H())
                    )
                );

                if (HVmfs.set(movingPhasei))
                {
                    phiHs[movingPhasei] += HVmfs[movingPhasei];
                }
            }

            phiHbyADs = invADVfs & phiHs;
        }

        forAll(movingPhases, movingPhasei)
        {
            const phaseModel& phase = movingPhases[movingPhasei];

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
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
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

        // Construct pressure "diffusivity"
        surfaceScalarField rAf
        (
            IOobject
            (
                "rAf",
                runTime.name(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimensionSet(-1, 3, 1, 0, 0), 0)
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
                    const label phasei = phase.index();

                    phase.phiRef() =
                        phiHbyADs[movingPhasei]
                      + alphaByADfs[movingPhasei]*mSfGradp;

                    // Set the phase dilatation rate
                    phase.divU(-pEqnComps[phasei] & p_rgh);
                }

                forAll(movingPhases, movingPhasei)
                {
                    phaseModel& phase = movingPhases_[movingPhasei];

                    MRF.makeRelative(phase.phiRef());
                    fvc::makeRelative(phase.phiRef(), phase.U());

                    phase.URef() = fvc::reconstruct
                    (
                        fvc::absolute(MRF.absolute(phase.phi()), phase.U())
                    );

                    phase.URef().correctBoundaryConditions();
                    phase.correctUf();
                    fvConstraints().constrain(phase.URef());
                }
            }
        }

        // Update and limit the static pressure
        p = p_rgh + rho*buoyancy.gh + buoyancy.pRef;
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
        p_rgh.correctBoundaryConditions();
    }

    UEqns.clear();
}


// ************************************************************************* //
