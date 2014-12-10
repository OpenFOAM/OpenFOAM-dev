{
    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1)
      - fvc::ddt(alpha1)
      - fvm::laplacian
        (
            volScalarField("Dab", Dab + alphatab*turbulence->nut()),
            alpha1
        )
    );

    alpha1Eqn.solve();

    alpha2 = 1.0 - alpha1;
    rhoPhi += alpha1Eqn.flux()*(rho1 - rho2);
}

rho = alpha1*rho1 + alpha2*rho2;
