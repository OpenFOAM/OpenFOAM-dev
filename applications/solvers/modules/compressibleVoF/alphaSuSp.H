tmp<volScalarField::Internal> Sp
(
    volScalarField::Internal::New
    (
        "Sp",
        mesh,
        dimensionedScalar(dgdt.dimensions(), 0)
    )
);

tmp<volScalarField::Internal> Su
(
    volScalarField::Internal::New
    (
        "Su",
        mesh,
        dimensionedScalar(dgdt.dimensions(), 0)
    )
);

if (fvModels().addsSupToField(alpha1.name()))
{
    // Phase change alpha1 source
    const fvScalarMatrix alphaSup(fvModels().source(alpha1));

    Su = alphaSup.Su();
    Sp = alphaSup.Sp();
}

volScalarField::Internal& SpRef = Sp.ref();
volScalarField::Internal& SuRef = Su.ref();

forAll(dgdt, celli)
{
    if (dgdt[celli] > 0.0)
    {
        SpRef[celli] -= dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
        SuRef[celli] += dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
    }
    else if (dgdt[celli] < 0.0)
    {
        SpRef[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
    }
}

tmp<volScalarField> divU
(
    mesh.moving()
  ? fvc::div(phiCN() + mesh.phi())
  : fvc::div(phiCN())
);
