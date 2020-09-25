{
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix eEqn
        (
            fvm::ddt(rho, e)
          + thermo.divq(e)
          ==
            fvOptions(rho, e)
        );

        eEqn.relax();

        fvOptions.constrain(eEqn);

        eEqn.solve();

        fvOptions.correct(e);
    }
}

thermo.correct();

Info<< "Min/max T:" << min(thermo.T()).value() << ' '
    << max(thermo.T()).value() << endl;
