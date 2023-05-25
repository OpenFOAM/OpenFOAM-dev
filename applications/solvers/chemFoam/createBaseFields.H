// write base thermo fields - not registered since will be re-read by
// thermo package

Info<< "Creating base fields for time " << runTime.name() << endl;
{
    volScalarField Ydefault
    (
        IOobject
        (
            "Ydefault",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimless, 1)
    );

    Ydefault.write();

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimPressure, p0)
    );

    p.write();

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimTemperature, T0)
    );

    T.write();
}
