{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he)
      + (
            mvConvection.valid()
          ? mvConvection->fvmDiv(phi, he)
          : fvm::div(phi, he)
        )
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      + thermophysicalTransport->divq(he)
     ==
        rho*(U&g)
      + clouds.Sh(he)
      + surfaceFilm.Sh()
      + combustion->Qdot()
      + fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
