if (ign.ignited())
{
    volScalarField& heau = thermo.heu();

    fvScalarMatrix heauEqn
    (
        betav*fvm::ddt(rho, heau) + mvConvection->fvmDiv(phi, heau)
      + (betav*fvc::ddt(rho, K) + fvc::div(phi, K))*rho/thermo.rhou()
      + (
            heau.name() == "eau"
          ? fvc::div(fvc::absolute(phi, rho, U), p/rho)
           *rho/thermo.rhou()
          : -betav*dpdt*rho/thermo.rhou()
        )
      - fvm::laplacian(Db, heau)

        // These terms cannot be used in partially-premixed combustion due to
        // the resultant inconsistency between ft and heau transport.
        // A possible solution would be to solve for ftu as well as ft.
        //- fvm::div(muEff*fvc::grad(b)/(b + 0.001), heau)
        //+ fvm::Sp(fvc::div(muEff*fvc::grad(b)/(b + 0.001)), heau)
     ==
        betav*fvModels.source(rho, heau)
    );

    fvConstraints.constrain(heauEqn);

    heauEqn.solve();

    fvConstraints.constrain(heau);
}
