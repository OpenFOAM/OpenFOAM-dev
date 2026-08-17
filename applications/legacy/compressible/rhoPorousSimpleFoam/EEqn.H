{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvi::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvi::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      + thermophysicalTransport->divq(he)
     ==
        fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();
}
