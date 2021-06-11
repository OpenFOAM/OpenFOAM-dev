{
    volScalarField& hea = thermo.he();

    fvScalarMatrix EaEqn
    (
        betav*fvm::ddt(rho, hea) + mvConvection->fvmDiv(phi, hea)
      + betav*fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            hea.name() == "ea"
          ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
          : -betav*dpdt
        )
      - fvm::laplacian(Db, hea)
      + betav*fvModels.source(rho, hea)
    );

    EaEqn.relax();

    fvConstraints.constrain(EaEqn);

    EaEqn.solve();

    fvConstraints.constrain(hea);

    thermo.correct();
}
