    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        betav*fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevTau(U)
     ==
        betav*rho*g
      + betav*fvOptions(rho, U)
    );

    fvOptions.constrain(UEqn);

    volSymmTensorField invA(inv(I*UEqn.A() + drag->Dcu()));

    if (pimple.momentumPredictor())
    {
        U = invA & (UEqn.H() - betav*fvc::grad(p));
        U.correctBoundaryConditions();
        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
