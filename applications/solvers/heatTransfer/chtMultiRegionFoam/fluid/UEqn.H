    // Solve the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turb.divDevRhoReff(U)
       ==
        fvOptions(rho, U)
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    if (momentumPredictor)
    {
        solve
        (
            UEqn()
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            ),
            mesh.solver(U.select(finalIter))
        );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
