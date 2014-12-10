if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        tmp<fvScalarMatrix> hEqn
        (
            fvm::ddt(betav*rho, h)
          - (
               thermo.isotropic()
             ? fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
             : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
            )
          ==
            fvOptions(rho, h)
        );

        hEqn().relax();

        fvOptions.constrain(hEqn());

        hEqn().solve(mesh.solver(h.select(finalIter)));

        fvOptions.correct(h);
    }
}

thermo.correct();

Info<< "Min/max T:" << min(thermo.T()) << ' ' << max(thermo.T()) << endl;

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
