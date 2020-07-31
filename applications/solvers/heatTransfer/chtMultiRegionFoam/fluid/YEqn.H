{
    reaction.correct();

    tmp<volScalarField> Yt(nullptr);

    if (Y.size())
    {
        Yt = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject("Yt", runTime.timeName(), mesh),
                mesh,
                dimensionedScalar(dimless, 0)
            )
        );
    }

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
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

            Yi.max(0.0);
            Yt.ref() += Yi;
        }
    }

    if (Y.size())
    {
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].max(0.0);
    }
}
