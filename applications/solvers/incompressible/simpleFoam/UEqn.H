    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevSigma(U)
     ==
        fvModels.source(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvConstraints.constrain(U);
    }
