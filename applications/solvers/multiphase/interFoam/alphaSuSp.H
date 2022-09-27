tmp<volScalarField> divU;
tmp<volScalarField::Internal> Su;
tmp<volScalarField::Internal> Sp;

if (fvModels.addsSupToField(alpha1.name()))
{
    // Phase change alpha1 source
    const fvScalarMatrix alphaSup(fvModels.source(alpha1));

    Su = alphaSup.Su();
    Sp = alphaSup.Sp();

    divU =
    (
        mesh.moving()
      ? fvc::div(phiCN() + mesh.phi())
      : fvc::div(phiCN())
    );
}
