if (thermo.containsSpecie("ft"))
{
    volScalarField& ft = thermo.Y("ft");

    solve
    (
        betav*fvm::ddt(rho, ft)
      + mvConvection->fvmDiv(phi, ft)
      - fvm::laplacian(Db, ft)
    );
}
