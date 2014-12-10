    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
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
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
