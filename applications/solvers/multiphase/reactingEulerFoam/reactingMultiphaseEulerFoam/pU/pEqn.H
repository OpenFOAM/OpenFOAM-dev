PtrList<surfaceScalarField> alphafs(phases.size());
PtrList<volScalarField> rAUs(phases.size());
PtrList<surfaceScalarField> alpharAUfs(phases.size());

// Diagonal coefficients
forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];
    const volScalarField& alpha = phase;

    alphafs.set(phasei, fvc::interpolate(alpha).ptr());
    alphafs[phasei].rename("pEqn" + alphafs[phasei].name());

    rAUs.set
    (
        phasei,
        new volScalarField
        (
            IOobject::groupName("rAU", phase.name()),
            1.0
           /(
               UEqns[phasei].A()
             + byDt(max(phase.residualAlpha() - alpha, scalar(0))*phase.rho())
            )
        )
    );

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

// --- Pressure corrector loop
while (pimple.correct())
{
    volScalarField rho("rho", fluid.rho());

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;

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

        forAll(phases, phasei)
        {
            phaseModel& phase = phases[phasei];
            const volScalarField& alpha = phase;

            // Correct fixed-flux BCs to be consistent with the velocity BCs
            MRF.correctBoundaryFlux(phase.U(), phase.phi());

            HbyAs.set
            (
                phasei,
                new volVectorField
                (
                    IOobject::groupName("HbyA", phase.name()),
                    phase.U()
                )
            );

            HbyAs[phasei] =
                rAUs[phasei]
               *(
                    UEqns[phasei].H()
                  + byDt
                    (
                        max(phase.residualAlpha() - alpha, scalar(0))
                       *phase.rho()
                    )
                   *phase.U().oldTime()
                );

            phiHbyAs.set
            (
                phasei,
                new surfaceScalarField
                (
                    IOobject::groupName("phiHbyA", phase.name()),
                    fvc::flux(HbyAs[phasei])
                  - phigFs[phasei]
                  - ddtCorrByAs[phasei]
                )
            );
        }
    }

    // Add explicit drag forces and fluxes if not doing partial elimination
    if (!partialElimination)
    {
        PtrList<volVectorField> KdUByAs(fluid.KdUByAs(rAUs));
        PtrList<surfaceScalarField> phiKdPhis(fluid.phiKdPhis(rAUs));

        forAll(phases, phasei)
        {
            HbyAs[phasei] -= KdUByAs[phasei];
            phiHbyAs[phasei] -= phiKdPhis[phasei];
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
        dimensionedScalar("phiHbyA", dimArea*dimVelocity, 0)
    );

    forAll(phases, phasei)
    {
        phiHbyA += alphafs[phasei]*phiHbyAs[phasei];
    }

    // Add explicit drag fluxes if doing partial elimination
    if (partialElimination)
    {
        PtrList<surfaceScalarField> phiKdPhis(fluid.phiKdPhis(rAUs));

        forAll(phases, phasei)
        {
            phiHbyA -= alphafs[phasei]*phiKdPhis[phasei];
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
        dimensionedScalar("rAUf", dimensionSet(-1, 3, 1, 0, 0), 0)
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
            phib += alphafs[phasei].boundaryField()*phase.phi().boundaryField();
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
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;
        volScalarField& rho = phase.thermo().rho();

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
                            phase.continuityError()
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

                pEqnComps[phasei].relax();
            }
            else
            {
                pEqnComps.set
                (
                    phasei,
                    (
                        (
                            phase.continuityError()
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
        else
        {
            pEqnComps.set
            (
                phasei,
                fvm::Su(-(fvOptions(alpha, rho)&rho)/rho, p_rgh).ptr()
            );
        }

        if (fluid.transfersMass())
        {
            if (pEqnComps.set(phasei))
            {
                pEqnComps[phasei] -= fluid.dmdt(phase)/rho;
            }
            else
            {
                pEqnComps.set
                (
                    phasei,
                    fvm::Su(-fluid.dmdt(phase)/rho, p_rgh)
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

            solve
            (
                pEqn,
                mesh.solver(p_rgh.select(pimple.finalInnerIter()))
            );
        }

        // Correct fluxes and velocities on last non-orthogonal iteration
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqnIncomp.flux();

            surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

            forAll(phases, phasei)
            {
                phaseModel& phase = phases[phasei];

                phase.phi() = phiHbyAs[phasei] + alpharAUfs[phasei]*mSfGradp;

                // Set the phase dilatation rates
                if (pEqnComps.set(phasei))
                {
                    phase.divU(-pEqnComps[phasei] & p_rgh);
                }
            }

            // Optionally relax pressure for velocity correction
            p_rgh.relax();

            mSfGradp = pEqnIncomp.flux()/rAUf;

            forAll(phases, phasei)
            {
                phaseModel& phase = phases[phasei];

                phase.U() =
                    HbyAs[phasei]
                  + fvc::reconstruct
                    (
                        alpharAUfs[phasei]*mSfGradp
                      - phigFs[phasei]
                    );
            }

            if (partialElimination)
            {
                fluid.partialElimination(rAUs);
            }

            forAll(phases, phasei)
            {
                phaseModel& phase = phases[phasei];

                phase.U().correctBoundaryConditions();
                fvOptions.correct(phase.U());
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
        phase.thermo().rho() += phase.thermo().psi()*(p_rgh - p_rgh_0);
    }

    // Correct p_rgh for consistency with p and the updated densities
    rho = fluid.rho();
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();
}
