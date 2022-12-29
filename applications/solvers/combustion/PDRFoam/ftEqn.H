if (composition.contains("ft"))
{
    volScalarField& ft = composition.Y("ft");

    solve
    (
        betav*fvm::ddt(rho, ft)
      + mvConvection->fvmDiv(phi, ft)
      - fvm::laplacian(Db, ft)
    );
}
