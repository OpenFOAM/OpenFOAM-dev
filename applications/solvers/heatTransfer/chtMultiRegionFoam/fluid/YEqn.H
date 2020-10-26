{
    reaction.correct();

    forAll(Y, i)
    {
        if (composition.solve(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi) + mvConvection->fvmDiv(phi, Yi)
              + thermophysicalTransport.divj(Yi)
             ==
                reaction.R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);
        }
    }

    composition.normalise();
}
