    // Relative momentum predictor

    tmp<fvVectorMatrix> UrelEqn
    (
        fvm::div(phi, Urel)
      + turbulence->divDevReff(Urel)
      + SRF->Su()
     ==
        fvOptions(Urel)
    );

    UrelEqn().relax();

    fvOptions.constrain(UrelEqn());

    solve(UrelEqn() == -fvc::grad(p));

    fvOptions.correct(Urel);
