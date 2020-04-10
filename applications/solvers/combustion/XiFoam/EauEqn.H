if (ign.ignited())
{
    volScalarField& heau = thermo.heu();

    fvScalarMatrix heauEqn
    (
        fvm::ddt(rho, heau) + mvConvection->fvmDiv(phi, heau)
      + (fvc::ddt(rho, K) + fvc::div(phi, K))*rho/thermo.rhou()
      + (
            heau.name() == "eau"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )*rho/thermo.rhou()
          : -dpdt*rho/thermo.rhou()
        )
      - fvm::laplacian(thermophysicalTransport->alphaEff(), heau)

        // These terms cannot be used in partially-premixed combustion due to
        // the resultant inconsistency between ft and heau transport.
        // A possible solution would be to solve for ftu as well as ft.
        //- fvm::div(muEff*fvc::grad(b)/(b + 0.001), heau)
        //+ fvm::Sp(fvc::div(muEff*fvc::grad(b)/(b + 0.001)), heau)

     ==
        fvOptions(rho, heau)
    );

    fvOptions.constrain(heauEqn);

    heauEqn.solve();

    fvOptions.correct(heau);
}
