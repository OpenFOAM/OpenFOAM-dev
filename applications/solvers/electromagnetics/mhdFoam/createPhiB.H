    IOobject phiBHeader
    (
        "phiB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );


    surfaceScalarField* phiBPtr;

    if (phiBHeader.headerOk())
    {
        Info<< "Reading field phiB\n" << endl;

        phiBPtr = new surfaceScalarField
        (
            IOobject
            (
                "phiB",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
    }
    else
    {
        Info<< "Calculating face flux\n" << endl;

        phiBPtr = new surfaceScalarField
        (
            IOobject
            (
                "phiB",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (fvc::interpolate(B) & mesh.Sf())
        );
    }

    surfaceScalarField& phiB = *phiBPtr;
