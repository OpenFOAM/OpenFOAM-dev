    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevTau(U)
     ==
        fvModels.source(rho, U)
    );

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvConstraints.constrain(U);
        K = 0.5*magSqr(U);
    }
