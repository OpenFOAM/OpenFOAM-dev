{
    bool closedVolume = p_rgh.needReference();
    dimensionedScalar compressibility = fvc::domainIntegrate(psi);
    bool compressible = (compressibility.value() > SMALL);

    rho = thermo.rho();

    volScalarField rAU("rAU", 1.0/UEqn().A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();

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

    {
        fvScalarMatrix p_rghDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
          + fvc::div(phiHbyA)
        );

        // Thermodynamic density needs to be updated by psi*d(p) after the
        // pressure solution - done in 2 parts. Part 1:
        thermo.rho() -= psi*p_rgh;

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix p_rghEqn
            (
                p_rghDDtEqn
              - fvm::laplacian(rhorAUf, p_rgh)
            );

            p_rghEqn.solve
            (
                mesh.solver
                (
                    p_rgh.select
                    (
                        (
                           oCorr == nOuterCorr-1
                        && corr == nCorr-1
                        && nonOrth == nNonOrthCorr
                        )
                    )
                )
            );

            if (nonOrth == nNonOrthCorr)
            {
                phi = phiHbyA + p_rghEqn.flux();
                U = HbyA
                  + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
                U.correctBoundaryConditions();
                fvOptions.correct(U);
                K = 0.5*magSqr(U);
            }
        }

        // Second part of thermodynamic density update
        thermo.rho() += psi*p_rgh;
    }

    p = p_rgh + rho*gh;

    // Update pressure time derivative if needed
    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);
    }

    // Solve continuity
    #include "rhoEqn.H"

    // Update continuity errors
    #include "compressibleContinuityErrors.H"

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume && compressible)
    {
        p += (initialMass - fvc::domainIntegrate(thermo.rho()))
            /compressibility;
        rho = thermo.rho();
        p_rgh = p - rho*gh;
    }
}
