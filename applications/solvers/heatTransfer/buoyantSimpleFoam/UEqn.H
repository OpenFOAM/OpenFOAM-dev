    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    if (simple.momentumPredictor())
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
            )
        );

        fvOptions.correct(U);
    }
