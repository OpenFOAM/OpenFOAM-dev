tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    combustion->correct();

    forAll(Y, i)
    {
        if (composition.solve(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YEqn
            (
                mvConvection->fvmDiv(phi, Yi)
              + thermophysicalTransport->divj(Yi)
             ==
                clouds.SYi(i, Yi)
              + combustion->R(Yi)
              + fvModels.source(rho, Yi)
            );

            YEqn.relax();

            fvConstraints.constrain(YEqn);

            YEqn.solve("Yi");

            fvConstraints.constrain(Yi);
        }
    }

    composition.normalise();
}
