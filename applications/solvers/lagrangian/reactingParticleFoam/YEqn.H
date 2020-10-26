tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);

if (Y.size())
{
    mvConvection = tmp<fv::convectionScheme<scalar>>
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );
}

{
    combustion->correct();

    forAll(Y, i)
    {
        if (composition.solve(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              + thermophysicalTransport->divj(Yi)
              ==
                clouds.SYi(i, Yi)
              + fvOptions(rho, Yi)
              + combustion->R(Yi)
              + surfaceFilm.Srho(i)
            );

            YEqn.relax();

            fvOptions.constrain(YEqn);

            YEqn.solve("Yi");

            fvOptions.correct(Yi);
        }
    }

    composition.normalise();
}
