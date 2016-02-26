    // Relative momentum predictor

    tmp<fvVectorMatrix> tUrelEqn
    (
        fvm::div(phi, Urel)
      + turbulence->divDevReff(Urel)
      + SRF->Su()
     ==
        fvOptions(Urel)
    );
    fvVectorMatrix& UrelEqn = tUrelEqn.ref();

    UrelEqn.relax();

    fvOptions.constrain(UrelEqn);

    solve(UrelEqn == -fvc::grad(p));

    fvOptions.correct(Urel);
