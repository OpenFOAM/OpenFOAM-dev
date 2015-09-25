PtrList<surfaceScalarField> alphafs(phases.size());
PtrList<volScalarField> rAUs(phases.size());
PtrList<surfaceScalarField> alpharAUfs(phases.size());

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
             + max(phase.residualAlpha() - alpha, scalar(0))
              *phase.rho()/runTime.deltaT()
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

// Lift, wall-lubrication and turbulent diffusion fluxes
PtrList<surfaceScalarField> phiFs(phases.size());
{
    autoPtr<PtrList<volVectorField> > Fs = fluid.Fs();

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        if (Fs().set(phasei))
        {
            phiFs.set
            (
                phasei,
                new surfaceScalarField
                (
                    IOobject::groupName("phiF", phase.name()),
                    (fvc::interpolate(rAUs[phasei]*Fs()[phasei]) & mesh.Sf())
                )
            );
        }
    }
}
{
    autoPtr<PtrList<surfaceScalarField> > phiDs = fluid.phiDs(rAUs);

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        if (phiDs().set(phasei))
        {
            if (phiFs.set(phasei))
            {
                phiFs[phasei] += phiDs()[phasei];
            }
            else
            {
                phiFs.set
                (
                    phasei,
                    new surfaceScalarField
                    (
                        IOobject::groupName("phiF", phase.name()),
                        phiDs()[phasei]
                    )
                );
            }
        }
    }
}

// --- Pressure corrector loop
while (pimple.correct())
{
    // Update continuity errors due to temperature changes
    fluid.correct();

    PtrList<volVectorField> HbyAs(phases.size());

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
              + max(phase.residualAlpha() - alpha, scalar(0))
               *phase.rho()*phase.U().oldTime()/runTime.deltaT()
            );
    }

    // Mean density for buoyancy force and p_rgh -> p
    volScalarField rho("rho", fluid.rho());

    surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    PtrList<surfaceScalarField> phigFs(phases.size());
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

    PtrList<surfaceScalarField> phiHbyAs(phases.size());

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
        phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;

        // ddtPhiCorr filter -- only apply in pure(ish) phases
        surfaceScalarField alphafBar
        (
            fvc::interpolate(fvc::average(alphafs[phasei]))
        );
        surfaceScalarField phiCorrCoeff(pos(alphafBar - 0.99));

        forAll(mesh.boundary(), patchi)
        {
            // Set ddtPhiCorr to 0 on non-coupled boundaries
            if
            (
                !mesh.boundary()[patchi].coupled()
             || isA<cyclicAMIFvPatch>(mesh.boundary()[patchi])
            )
            {
                phiCorrCoeff.boundaryField()[patchi] = 0;
            }
        }

        phiHbyAs.set
        (
            phasei,
            new surfaceScalarField
            (
                IOobject::groupName("phiHbyA", phase.name()),
                (fvc::interpolate(HbyAs[phasei]) & mesh.Sf())
              + phiCorrCoeff
               *fvc::interpolate
                (
                    alpha.oldTime()*phase.rho()().oldTime()*rAUs[phasei]
                )
               *(
                   MRF.absolute(phase.phi().oldTime())
                 - (fvc::interpolate(phase.U().oldTime()) & mesh.Sf())
                )/runTime.deltaT()
              - phigFs[phasei]
            )
        );

        forAllConstIter
        (
            phaseSystem::KdTable,
            fluid.Kds(),
            KdIter
        )
        {
            const volScalarField& K(*KdIter());

            const phasePair& pair(fluid.phasePairs()[KdIter.key()]);

            const phaseModel* phase1 = &pair.phase1();
            const phaseModel* phase2 = &pair.phase2();

            forAllConstIter(phasePair, pair, iter)
            {
                if (phase1 == &phase)
                {
                    phiHbyAs[phasei] +=
                        fvc::interpolate(rAUs[phasei]*K)
                       *MRF.absolute(phase2->phi());

                    HbyAs[phasei] += rAUs[phasei]*K*phase2->U();
                }

                Swap(phase1, phase2);
            }
        }

        phiHbyA += alphafs[phasei]*phiHbyAs[phasei];
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
        surfaceScalarField::GeometricBoundaryField phib(phi.boundaryField());
        phib = 0;
        forAll(phases, phasei)
        {
            phaseModel& phase = phases[phasei];
            phib += alphafs[phasei].boundaryField()*phase.phi().boundaryField();
        }

        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p_rgh.boundaryField(),
            (
                phiHbyA.boundaryField() - phib
            )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
        );
    }

    PtrList<fvScalarMatrix> pEqnComps(phases.size());
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        if (phase.compressible())
        {
            const volScalarField& alpha = phase;
            const volScalarField& rho = phase.rho();

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
                      + (alpha/rho)*correction
                        (
                            phase.thermo().psi()*fvm::ddt(p_rgh)
                          + fvm::div(phid, p_rgh)
                          - fvm::Sp(fvc::div(phid), p_rgh)
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

            if (fluid.transfersMass(phase))
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
                if (phase.compressible())
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
        phase.rho()() += phase.thermo().psi()*(p_rgh - p_rgh_0);
    }

    // Correct p_rgh for consistency with p and the updated densities
    rho = fluid.rho();
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();
}
