{
    volScalarField& hea = thermo.he();

    solve
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
    );

    thermo.correct();
}
