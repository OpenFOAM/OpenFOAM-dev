{
    volScalarField& he1 = thermo1.he();
    volScalarField& he2 = thermo2.he();

    volScalarField Cpv1("Cpv1", thermo1.Cpv());
    volScalarField Cpv2("Cpv2", thermo2.Cpv());

    volScalarField Kh(fluid.Kh());

    fvScalarMatrix E1Eqn
    (
        fvm::ddt(alpha1, rho1, he1) + fvm::div(alphaRhoPhi1, he1)
      - fvm::Sp(contErr1, he1)

      + fvc::ddt(alpha1, rho1, K1) + fvc::div(alphaRhoPhi1, K1)
      - contErr1*K1
      + (
            he1.name() == thermo1.phasePropertyName("e")
          ? fvc::div(fvc::absolute(alphaPhi1, alpha1, U1), p)
          + p*fvc::ddt(alpha1)
          : -alpha1*dpdt
        )

      - fvm::laplacian
        (
            fvc::interpolate(alpha1)
           *fvc::interpolate(thermo1.alphaEff(phase1.turbulence().mut())),
            he1
        )
    );

    E1Eqn.relax();

    E1Eqn -=
    (
        Kh*(thermo2.T() - thermo1.T())
      + Kh*he1/Cpv1
      - fvm::Sp(Kh/Cpv1, he1)
      + alpha1*rho1*(U1&g)
      + fvOptions(alpha1, rho1, he1)
    );

    fvScalarMatrix E2Eqn
    (
        fvm::ddt(alpha2, rho2, he2) + fvm::div(alphaRhoPhi2, he2)
      - fvm::Sp(contErr2, he2)

      + fvc::ddt(alpha2, rho2, K2) + fvc::div(alphaRhoPhi2, K2)
      - contErr2*K2
      + (
            he2.name() == thermo2.phasePropertyName("e")
          ? fvc::div(fvc::absolute(alphaPhi2, alpha2, U2), p)
          + p*fvc::ddt(alpha1)
          : -alpha2*dpdt
        )

      - fvm::laplacian
        (
            fvc::interpolate(alpha2)
           *fvc::interpolate(thermo2.alphaEff(phase2.turbulence().mut())),
            he2
        )
    );

    E2Eqn.relax();

    E2Eqn -=
    (
        Kh*(thermo1.T() - thermo2.T())
      + Kh*he2/Cpv2
      - fvm::Sp(Kh/Cpv2, he2)
      + alpha2*rho2*(U2&g)
      + fvOptions(alpha2, rho2, he2)
    );

    fvOptions.constrain(E1Eqn);
    E1Eqn.solve();
    fvOptions.correct(he1);

    fvOptions.constrain(E2Eqn);
    E2Eqn.solve();
    fvOptions.correct(he2);

    thermo1.correct();
    Info<< "min " << thermo1.T().name()
        << " " << min(thermo1.T()).value() << endl;

    thermo2.correct();
    Info<< "min " << thermo2.T().name()
        << " " << min(thermo2.T()).value() << endl;
}
