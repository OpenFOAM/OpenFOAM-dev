{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    PtrList<fvScalarMatrix> p_rghEqnComps(mixture.phases().size());

    forAll(mixture.phases(), phasei)
    {
        const compressiblePhase& phase = mixture.phases()[phasei];
        const rhoThermo& thermo = phase.thermo();
        const volScalarField& rho = thermo.rho()();

        p_rghEqnComps.set
        (
            phasei,
            (
                fvc::ddt(rho) + thermo.psi()*correction(fvm::ddt(p_rgh))
              + fvc::div(phi, rho) - fvc::Sp(fvc::div(phi), rho)
              - (fvModels.source(phase, rho)&rho)
            ).ptr()
        );
    }

    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        tmp<fvScalarMatrix> p_rghEqnComp;

        forAll(mixture.phases(), phasei)
        {
            const compressiblePhase& phase = mixture.phases()[phasei];

            tmp<fvScalarMatrix> p_rghEqnCompi
            (
                (max(phase, scalar(0))/phase.thermo().rho())
               *p_rghEqnComps[phasei]
            );

            if (phasei == 0)
            {
                p_rghEqnComp = p_rghEqnCompi;
            }
            else
            {
                p_rghEqnComp.ref() += p_rghEqnCompi;
            }
        }

        solve(p_rghEqnComp + p_rghEqnIncomp);

        if (pimple.finalNonOrthogonalIter())
        {
            p = max(p_rgh + mixture.rho()*gh, pMin);
            p_rgh = p - rho*gh;

            forAll(mixture.phases(), phasei)
            {
                compressiblePhase& phase = mixture.phases()[phasei];

                phase.dgdt() =
                    pos0(phase)
                  *(p_rghEqnComps[phasei] & p_rgh)/phase.thermo().rho();
            }

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
            fvConstraints.constrain(U);
        }
    }

    // Update densities from change in p_rgh
    mixture.correctRho(p_rgh - p_rgh_0);
    mixture.correct();

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();

    K = 0.5*magSqr(U);

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;
}
