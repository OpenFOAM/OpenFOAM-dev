{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix hEqn
        (
            (
                thermo.isotropic()
                ? -fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
                : -fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
            )
          ==
            fvOptions(rho, h)
        );

        hEqn.relax();

        fvOptions.constrain(hEqn);

        hEqn.solve();

        fvOptions.correct(h);
    }
}

thermo.correct();

Info<< "Min/max T:" << min(thermo.T()) << ' ' << max(thermo.T()) << endl;
