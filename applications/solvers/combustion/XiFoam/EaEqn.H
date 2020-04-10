{
    volScalarField& hea = thermo.he();

    fvScalarMatrix EaEqn
    (
        fvm::ddt(rho, hea) + mvConvection->fvmDiv(phi, hea)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            hea.name() == "ea"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(thermophysicalTransport->alphaEff(), hea)
      + fvOptions(rho, hea)
    );

    EaEqn.relax();

    fvOptions.constrain(EaEqn);

    EaEqn.solve();

    fvOptions.correct(hea);

    thermo.correct();
}
