tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,ft_b_ha_hau)")
    )
);

if (composition.contains("ft"))
{
    volScalarField& ft = composition.Y("ft");

    fvScalarMatrix ftEqn
    (
        fvm::ddt(rho, ft)
      + mvConvection->fvmDiv(phi, ft)
      - fvm::laplacian(turbulence->alphaEff(), ft)
     ==
        fvOptions(rho, ft)
    );

    fvOptions.constrain(ftEqn);

    ftEqn.solve();

    fvOptions.correct(ft);
}
