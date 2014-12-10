{
    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p;

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            (fvc::interpolate(rho*HbyA) & mesh.Sf())
          + rhorAUf*fvc::ddtCorr(rho, U, phi)
        )
      + phig
    );

    fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
           *rho.boundaryField()
        )/(mesh.magSf().boundaryField()*rhorAUf.boundaryField())
    );

    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
      + fvc::div(phiHbyA)
      ==
        fvOptions(psi, p_rgh, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            p_rghDDtEqn
          - fvm::laplacian(rhorAUf, p_rgh)
        );

        fvOptions.constrain(p_rghEqn);

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi = phiHbyA + p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
            K = 0.5*magSqr(U);
        }
    }

    p = p_rgh + rho*gh;

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);
    }

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"
}
