// write base thermo fields - not registered since will be re-read by
// thermo package

Info<< "Creating base fields for time " << runTime.timeName() << endl;
{
    volScalarField Ydefault
    (
        IOobject
        (
            "Ydefault",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("Ydefault", dimless, 1)
    );

    Ydefault.write();

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("p", dimPressure, p0)
    );

    p.write();

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, T0)
    );

    T.write();
}

