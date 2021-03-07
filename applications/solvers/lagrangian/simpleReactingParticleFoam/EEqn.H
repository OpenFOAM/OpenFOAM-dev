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
