{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        mvConvection->fvmDiv(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      + thermophysicalTransport->divq(he)
     ==
        rho*(U&g)
      + clouds.Sh(he)
      + radiation->Sh(thermo, he)
      + combustion->Qdot()
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);
    thermo.correct();
    radiation->correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
