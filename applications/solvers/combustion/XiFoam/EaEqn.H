{
    volScalarField& hea = thermo.he();

    fvScalarMatrix EaEqn
    (
        fvm::ddt(rho, hea) + mvConvection->fvmDiv(phi, hea)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            hea.name() == "ea"
          ? fvc::div(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
      + thermophysicalTransport.divq(hea)
      + fvModels.source(rho, hea)
    );

    EaEqn.relax();

    fvConstraints.constrain(EaEqn);

    EaEqn.solve();

    fvConstraints.constrain(hea);

    thermo.correct();
}
