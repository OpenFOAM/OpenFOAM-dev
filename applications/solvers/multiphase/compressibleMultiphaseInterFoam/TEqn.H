{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(mixture.rhoPhi(), T)
      - fvm::laplacian(mixture.alphaEff(rho*turbulence->nut()), T)
      + (
            fvc::div(fvc::absolute(phi, U), p)
          + fvc::ddt(rho, K) + fvc::div(mixture.rhoPhi(), K)
        )*mixture.rCv()
    );

    TEqn.relax();
    TEqn.solve();

    mixture.correctThermo();
    mixture.correct();
}
