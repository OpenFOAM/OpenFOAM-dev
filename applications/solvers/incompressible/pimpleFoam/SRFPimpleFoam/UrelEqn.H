    // Relative momentum predictor
    tmp<fvVectorMatrix> tUrelEqn
    (
        fvm::ddt(Urel)
      + fvm::div(phi, Urel)
      + turbulence->divDevSigma(Urel)
      + SRF->Su()
     ==
        fvModels.source(Urel)
    );
    fvVectorMatrix& UrelEqn = tUrelEqn.ref();

    UrelEqn.relax();

    fvConstraints.constrain(UrelEqn);

    solve(UrelEqn == -fvc::grad(p));

    fvConstraints.constrain(Urel);
