{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvm::div(phi, rho)
    );

    rhoEqn.solve();

    rhoPhi = rhoEqn.flux();

    Info<< "max-min rho: " << max(rho).value()
        << " " << min(rho).value() << endl;

    rho == max(rho, rhoMin);
}
