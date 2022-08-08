    // Relative momentum predictor

    tmp<fvVectorMatrix> tUrelEqn
    (
        fvm::div(phi, Urel)
      + turbulence->divDevSigma(Urel)
      + SRF->Su()
     ==
        fvModels.source(Urel)
    );
    fvVectorMatrix& UrelEqn = tUrelEqn.ref();

    UrelEqn.relax();

    fvConstraints.constrain(UrelEqn);

    if (simple.momentumPredictor())
    {
        solve(UrelEqn == -fvc::grad(p));

        fvConstraints.constrain(Urel);
    }
