volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Sp", dgdt.dimensions(), 0)
);

volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Su", dgdt.dimensions(), 0)
);

forAll(dgdt, celli)
{
    if (dgdt[celli] > 0.0 && alpha1[celli] > 0.0)
    {
        Sp[celli] -= dgdt[celli]*alpha1[celli];
        Su[celli] += dgdt[celli]*alpha1[celli];
    }
    else if (dgdt[celli] < 0.0 && alpha1[celli] < 1.0)
    {
        Sp[celli] += dgdt[celli]*(1.0 - alpha1[celli]);
    }
}

volScalarField::Internal divU
(
    mesh.moving()
  ? fvc::div(phiCN() + mesh.phi())
  : fvc::div(phiCN())
);
