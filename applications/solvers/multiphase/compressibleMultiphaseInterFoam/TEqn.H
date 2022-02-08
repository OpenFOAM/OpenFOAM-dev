{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T) + fvm::div(mixture.rhoPhi(), T) - fvm::Sp(contErr, T)
      - fvm::laplacian(mixture.alphaEff(rho*turbulence->nut()), T)
      + (
            fvc::div(fvc::absolute(phi, U), p)()() // - contErr/rho*p
          + (fvc::ddt(rho, K) + fvc::div(mixture.rhoPhi(), K))()()
          - (U()&(fvModels.source(rho, U)&U)()) - contErr*K
        )*mixture.rCv()()
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
