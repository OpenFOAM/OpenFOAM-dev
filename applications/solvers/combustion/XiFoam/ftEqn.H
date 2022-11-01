tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.schemes().div("div(phi,ft_b_ha_hau)")
    )
);

if (composition.contains("ft"))
{
    volScalarField& ft = composition.Y("ft");

    fvScalarMatrix ftEqn
    (
        fvm::ddt(rho, ft)
      + mvConvection->fvmDiv(phi, ft)
      - fvm::laplacian(thermophysicalTransport.DEff(ft), ft)
     ==
        fvModels.source(rho, ft)
    );

    fvConstraints.constrain(ftEqn);

    ftEqn.solve();

    fvConstraints.constrain(ft);
}
