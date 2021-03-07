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

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              + thermophysicalTransport->divj(Yi)
              ==
                parcels.SYi(i, Yi)
              + surfaceFilm.Srho(i)
              + combustion->R(Yi)
              + fvModels.source(rho, Yi)
            );

            YiEqn.relax();

            fvConstraints.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvConstraints.constrain(Yi);
        }
    }

    composition.normalise();

    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      + thermophysicalTransport->divq(he)
     ==
        combustion->Qdot()
      + parcels.Sh(he)
      + surfaceFilm.Sh()
      + fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
