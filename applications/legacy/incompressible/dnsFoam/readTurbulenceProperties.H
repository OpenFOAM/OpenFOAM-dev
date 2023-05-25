    Info<< "Reading momentumTransport\n" << endl;

    IOdictionary momentumTransport
    (
        IOobject
        (
            "momentumTransport",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    volVectorField force
    (
        U/dimensionedScalar("dt", dimTime, runTime.deltaTValue())
    );

    Kmesh K(mesh);
    UOprocess forceGen(K, runTime.deltaTValue(), momentumTransport);
