    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        betav*fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U)
     ==
        betav*rho*g
    );

    volSymmTensorField invA(inv(I*UEqn.A() + drag->Dcu()));

    if (pimple.momentumPredictor())
    {
        U = invA & (UEqn.H() - betav*fvc::grad(p));
        U.correctBoundaryConditions();
        K = 0.5*magSqr(U);
    }
