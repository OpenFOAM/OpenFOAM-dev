{
    volScalarField& hea = thermo.he();

    fvScalarMatrix EaEqn
    (
        betav*fvm::ddt(rho, hea) + mvConvection->fvmDiv(phi, hea)
      + betav()*fvi::ddt(rho, K) + fvi::div(phi, K)
      + (
            hea.name() == "ea"
          ? fvi::div(fvc::absolute(phi, rho, U), p/rho)
          : -betav()*dpdt
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
