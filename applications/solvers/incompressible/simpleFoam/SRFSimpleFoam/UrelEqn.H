    // Relative momentum predictor

    tmp<fvVectorMatrix> tUrelEqn
    (
        fvm::div(phi, Urel)
      + turbulence->divDevSigma(Urel)
      + SRF->Su()
     ==
        fvOptions(Urel)
    );
    fvVectorMatrix& UrelEqn = tUrelEqn.ref();

    UrelEqn.relax();

    fvOptions.constrain(UrelEqn);

    if (simple.momentumPredictor())
    {
        solve(UrelEqn == -fvc::grad(p));

        fvOptions.correct(Urel);
    }
