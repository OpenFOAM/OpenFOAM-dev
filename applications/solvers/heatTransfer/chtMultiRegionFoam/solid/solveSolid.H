{
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix eEqn
        (
            fvm::ddt(rho, e)
          + thermo.divq(e)
          ==
            fvModels.source(rho, e)
        );

        eEqn.relax();

        fvConstraints.constrain(eEqn);

        eEqn.solve();

        fvConstraints.constrain(e);
    }
}

thermo.correct();

Info<< "Min/max T:" << min(thermo.T()).value() << ' '
    << max(thermo.T()).value() << endl;
