{
    volScalarField& he1 = thermo1.he();
    volScalarField& he2 = thermo2.he();

    volScalarField Cpv1("Cpv1", thermo1.Cpv());
    volScalarField Cpv2("Cpv2", thermo2.Cpv());

    volScalarField heatTransferCoeff(fluid.heatTransferCoeff());

    fvScalarMatrix he1Eqn
    (
        fvm::ddt(alpha1, rho1, he1) + fvm::div(alphaRhoPhi1, he1)
      - fvm::Sp(contErr1, he1)

      + fvc::ddt(alpha1, rho1, K1) + fvc::div(alphaRhoPhi1, K1)
      - contErr1*K1
      + (
            he1.name() == thermo1.phasePropertyName("e")
          ? fvc::ddt(alpha1)*p + fvc::div(alphaPhi1, p)
          : -alpha1*dpdt
        )

      - fvm::laplacian
        (
            fvc::interpolate(alpha1)
           *fvc::interpolate(thermo1.alphaEff(phase1.turbulence().mut())),
            he1
        )
    );

    he1Eqn.relax();

    he1Eqn -=
    (
        heatTransferCoeff*(thermo2.T() - thermo1.T())
      + heatTransferCoeff*he1/Cpv1
      - fvm::Sp(heatTransferCoeff/Cpv1, he1)
      + fvOptions(alpha1, rho1, he1)
    );

    fvScalarMatrix he2Eqn
    (
        fvm::ddt(alpha2, rho2, he2) + fvm::div(alphaRhoPhi2, he2)
      - fvm::Sp(contErr2, he2)

      + fvc::ddt(alpha2, rho2, K2) + fvc::div(alphaRhoPhi2, K2)
      - contErr2*K2
      + (
            he2.name() == thermo2.phasePropertyName("e")
          ? fvc::ddt(alpha2)*p + fvc::div(alphaPhi2, p)
          : -alpha2*dpdt
        )

      - fvm::laplacian
        (
            fvc::interpolate(alpha2)
           *fvc::interpolate(thermo2.alphaEff(phase2.turbulence().mut())),
            he2
        )
    );

    he2Eqn.relax();

    he2Eqn -=
    (
        heatTransferCoeff*(thermo1.T() - thermo2.T())
      + heatTransferCoeff*he2/Cpv2
      - fvm::Sp(heatTransferCoeff/Cpv2, he2)
      + fvOptions(alpha2, rho2, he2)
    );

    fvOptions.constrain(he1Eqn);
    he1Eqn.solve();

    fvOptions.constrain(he2Eqn);
    he2Eqn.solve();

    thermo1.correct();
    Info<< "min " << thermo1.T().name()
        << " " << min(thermo1.T()).value() << endl;

    thermo2.correct();
    Info<< "min " << thermo2.T().name()
        << " " << min(thermo2.T()).value() << endl;
}
