{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T) + fvm::div(rhoPhi, T) - fvm::Sp(contErr, T)
      - fvm::laplacian(turbulence.alphaEff(), T)
      + (
            fvc::div(fvc::absolute(phi, U), p)()() // - contErr/rho*p
          + (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))()()
          - (U()&(fvModels.source(rho, U)&U)()) - contErr*K
        )
       *(
           alpha1()/mixture.thermo1().Cv()()
         + alpha2()/mixture.thermo2().Cv()()
        )
     ==
        fvModels.source(rho, T)
    );

    TEqn.relax();

    fvConstraints.constrain(TEqn);

    TEqn.solve();

    fvConstraints.constrain(T);

    mixture.correctThermo();
    mixture.correct();
}
