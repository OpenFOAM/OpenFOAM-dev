{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvm::laplacian(mixture.alphaEff(turbulence->mut()), T)
      + (
            fvc::div(fvc::absolute(phi, U), p)
          + fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
        )
       *(
           alpha1/mixture.thermo1().Cv()
         + alpha2/mixture.thermo2().Cv()
        )
    );

    TEqn.relax();
    TEqn.solve();

    mixture.correct();

    Info<< "min(T) " << min(T).value() << endl;
}
