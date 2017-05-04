    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turb.divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (momentumPredictor)
    {
        solvPerfU =
            solve
            (
                UEqn
              ==
                fvc::reconstruct
                (
                    (
                      - ghf*fvc::snGrad(rho)
                      - fvc::snGrad(p_rgh)
                    )*mesh.magSf()
                ),
                mesh.solver(U.select(finalIter))
            );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
