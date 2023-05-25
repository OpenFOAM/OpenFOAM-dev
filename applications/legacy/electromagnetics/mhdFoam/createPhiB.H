typeIOobject<surfaceScalarField> phiBHeader
(
    "phiB",
    runTime.name(),
    mesh,
    IOobject::NO_READ
);


surfaceScalarField* phiBPtr = nullptr;

if (phiBHeader.headerOk())
{
    Info<< "Reading face flux ";

    phiBPtr = new surfaceScalarField
    (
        IOobject
        (
            "phiB",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
}
else
{
    Info<< "Calculating face flux ";

    phiBPtr = new surfaceScalarField
    (
        IOobject
        (
            "phiB",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(B)
    );
}

surfaceScalarField& phiB = *phiBPtr;

Info<< phiB.name() << nl << endl;
