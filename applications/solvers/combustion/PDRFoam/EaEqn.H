{
    volScalarField& hea = thermo.he();

    fvScalarMatrix EaEqn
    (
        betav*fvm::ddt(rho, hea) + mvConvection->fvmDiv(phi, hea)
      + betav*fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            hea.name() == "ea"
          ? fvc::div
            (
                phi/fvc::interpolate(rho),
                p,
                "div(phiv,p)"
            )
          : -betav*dpdt
        )
      - fvm::laplacian(Db, hea)
      + betav*fvOptions(rho, hea)
    );

    EaEqn.relax();

    fvOptions.constrain(EaEqn);

    EaEqn.solve();

    fvOptions.correct(hea);

    thermo.correct();
}
