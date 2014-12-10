    Info<< "Reading boxTurbDict\n" << endl;

    IOdictionary boxTurbDict
    (
        IOobject
        (
            "boxTurbDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar Ea(readScalar(boxTurbDict.lookup("Ea")));

    scalar k0(readScalar(boxTurbDict.lookup("k0")));
