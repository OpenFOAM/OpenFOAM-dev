// Face volume fractions
PtrList<surfaceScalarField> alphafs(phases.size());
PtrList<surfaceScalarField> alphaRho0fs(phases.size());
forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];
    const volScalarField& alpha = phase;

    alphafs.set(phasei, fvc::interpolate(alpha).ptr());
    alphafs[phasei].rename("pEqn" + alphafs[phasei].name());

    alphaRho0fs.set
    (
        phasei,
        (
            fvc::interpolate
            (
                max(alpha.oldTime(), phase.residualAlpha())
               *phase.rho()().oldTime()
            )
        ).ptr()
    );
}

// Diagonal coefficients
PtrList<surfaceScalarField> rAUfs(phases.size());
{
    PtrList<surfaceScalarField> AFfs(fluid.AFfs());

    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];

        rAUfs.set
        (
            phase.index(),
            new surfaceScalarField
            (
                IOobject::groupName("rAUf", phase.name()),
                1.0
               /(
                    byDt(alphaRho0fs[phase.index()])
                  + fvc::interpolate(UEqns[phase.index()].A())
                  + AFfs[phase.index()]
                )
            )
        );
    }
}
fluid.fillFields("rAUf", dimTime/dimDensity, rAUfs);

// Phase diagonal coefficients
PtrList<surfaceScalarField> alpharAUfs(phases.size());
forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];
    alpharAUfs.set
    (
        phase.index(),
        (
            max(alphafs[phase.index()], phase.residualAlpha())
           *rAUfs[phase.index()]
        ).ptr()
    );
}

// Explicit force fluxes
PtrList<surfaceScalarField> phiFfs(fluid.phiFfs(rAUfs));

// --- Pressure corrector loop
while (pimple.correct())
{
    volScalarField rho("rho", fluid.rho());

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];
        MRF.correctBoundaryFlux(phase.U(), phase.phiRef());
    }

    // Combined buoyancy and force fluxes
    PtrList<surfaceScalarField> phigFs(phases.size());
    {
        const surfaceScalarField ghSnGradRho
        (
            "ghSnGradRho",
            ghf*fvc::snGrad(rho)*mesh.magSf()
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
                     - (fvc::interpolate(phase.rho() - rho))*(g & mesh.Sf())
                     - fluid.surfaceTension(phase)*mesh.magSf()
                    )
                ).ptr()
            );

            if (phiFfs.set(phasei))
            {
                phigFs[phasei] += phiFfs[phasei];
            }
        }
    }

    // Predicted fluxes for each phase
    PtrList<surfaceScalarField> phiHbyAs(phases.size());
    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];

        phiHbyAs.set
        (
            phase.index(),
            new surfaceScalarField
            (
                IOobject::groupName("phiHbyA", phase.name()),
                rAUfs[phase.index()]
               *(
                    fvc::flux(UEqns[phase.index()].H())
                  + alphaRho0fs[phase.index()]
                   *byDt(MRF.absolute(phase.phi()().oldTime()))
                )
              - phigFs[phase.index()]
            )
        );
    }
    fluid.fillFields("phiHbyA", dimForce/dimDensity/dimVelocity, phiHbyAs);

    // Add explicit drag forces and fluxes if not doing partial elimination
    if (!partialElimination)
    {
        PtrList<surfaceScalarField> phiKdPhifs(fluid.phiKdPhifs(rAUfs));

        forAll(phases, phasei)
        {
            if (phiKdPhifs.set(phasei))
            {
                phiHbyAs[phasei] -= phiKdPhifs[phasei];
            }
        }
    }

    // Total predicted flux
    surfaceScalarField phiHbyA
    (
        IOobject
        (
            "phiHbyA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimArea*dimVelocity, 0)
    );

    forAll(phases, phasei)
    {
        phiHbyA += alphafs[phasei]*phiHbyAs[phasei];
    }

    // Add explicit drag fluxes if doing partial elimination
    if (partialElimination)
    {
        PtrList<surfaceScalarField> phiKdPhifs(fluid.phiKdPhifs(rAUfs));

        forAll(phases, phasei)
        {
            if (phiKdPhifs.set(phasei))
            {
                phiHbyA -= alphafs[phasei]*phiKdPhifs[phasei];
            }
        }
    }

    MRF.makeRelative(phiHbyA);

    // Construct pressure "diffusivity"
    surfaceScalarField rAUf
    (
        IOobject
        (
            "rAUf",
            runTime.timeName(),
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
        surfaceScalarField::Boundary phib(phi.boundaryField());
        phib = 0;
        forAll(phases, phasei)
        {
            phaseModel& phase = phases[phasei];
            phib +=
                alphafs[phasei].boundaryField()*phase.phi()().boundaryField();
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
    PtrList<fvScalarMatrix> pEqnComps(phases.size());
    PtrList<volScalarField> dmdts(fluid.dmdts());
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;
        volScalarField& rho = phase.thermoRef().rho();

        if (phase.compressible())
        {
            if (pimple.transonic())
            {
                surfaceScalarField phid
                (
                    IOobject::groupName("phid", phase.name()),
                    fvc::interpolate(phase.thermo().psi())*phase.phi()
                );

                pEqnComps.set
                (
                    phasei,
                    (
                        (
                            fvc::ddt(alpha, rho) + fvc::div(phase.alphaRhoPhi())
                          - fvc::Sp
                            (
                                fvc::ddt(alpha) + fvc::div(phase.alphaPhi()),
                                rho
                            )
                        )/rho
                      + correction
                        (
                            (alpha/rho)*
                            (
                                phase.thermo().psi()*fvm::ddt(p_rgh)
                              + fvm::div(phid, p_rgh)
                              - fvm::Sp(fvc::div(phid), p_rgh)
                            )
                        )
                    ).ptr()
                );

                deleteDemandDrivenData
                (
                    pEqnComps[phasei].faceFluxCorrectionPtr()
                );

                pEqnComps[phasei].relax();
            }
            else
            {
                pEqnComps.set
                (
                    phasei,
                    (
                        (
                            fvc::ddt(alpha, rho) + fvc::div(phase.alphaRhoPhi())
                          - fvc::Sp
                            (
                                (fvc::ddt(alpha) + fvc::div(phase.alphaPhi())),
                                rho
                            )
                        )/rho
                      + (alpha*phase.thermo().psi()/rho)
                       *correction(fvm::ddt(p_rgh))
                    ).ptr()
                );
            }
        }

        if (fvOptions.appliesToField(rho.name()))
        {
            tmp<fvScalarMatrix> optEqn = fvOptions(alpha, rho);
            if (pEqnComps.set(phasei))
            {
                pEqnComps[phasei] -= (optEqn&rho)/rho;
            }
            else
            {
                pEqnComps.set
                (
                    phasei,
                    fvm::Su(- (optEqn&rho)/rho, p_rgh).ptr()
                );
            }
        }

        if (dmdts.set(phasei))
        {
            if (pEqnComps.set(phasei))
            {
                pEqnComps[phasei] -= dmdts[phasei]/rho;
            }
            else
            {
                pEqnComps.set
                (
                    phasei,
                    fvm::Su(- dmdts[phasei]/rho, p_rgh)
                );
            }
        }
    }

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

        {
            fvScalarMatrix pEqn(pEqnIncomp);

            forAll(phases, phasei)
            {
                if (pEqnComps.set(phasei))
                {
                    pEqn += pEqnComps[phasei];
                }
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

                // Set the phase dilatation rates
                if (pEqnComps.set(phase.index()))
                {
                    phase.divU(-pEqnComps[phase.index()] & p_rgh);
                }
            }

            if (partialElimination)
            {
                fluid.partialEliminationf(rAUfs);
            }

            // Optionally relax pressure for velocity correction
            p_rgh.relax();

            mSfGradp = pEqnIncomp.flux()/rAUf;

            forAll(fluid.movingPhases(), movingPhasei)
            {
                phaseModel& phase = fluid.movingPhases()[movingPhasei];

                phase.URef() = fvc::reconstruct(MRF.absolute(phase.phi()));
                phase.URef().correctBoundaryConditions();
                fvOptions.correct(phase.URef());
            }
        }
    }

    // Update and limit the static pressure
    p = max(p_rgh + rho*gh, pMin);

    // Limit p_rgh
    p_rgh = p - rho*gh;

    // Update densities from change in p_rgh
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        phase.thermoRef().rho() += phase.thermo().psi()*(p_rgh - p_rgh_0);
    }

    // Correct p_rgh for consistency with p and the updated densities
    rho = fluid.rho();
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();
}
