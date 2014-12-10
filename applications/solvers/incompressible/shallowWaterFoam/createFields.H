Info<< "Reading field h\n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field h0 if present\n" << endl;
volScalarField h0
(
    IOobject
    (
        "h0",
        runTime.findInstance("polyMesh", "points"),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mesh,
    dimensionedScalar("h0", dimLength, 0.0)
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating field hU\n" << endl;
volVectorField hU
(
    IOobject
    (
        "hU",
        runTime.timeName(),
        mesh
    ),
    h*U,
    U.boundaryField().types()
);

Info<< "Creating field hTotal for post processing\n" << endl;
volScalarField hTotal
(
    IOobject
    (
        "hTotal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    h+h0
);
hTotal.write();

#   include "createPhi.H"

Info<< "Creating Coriolis Force" << endl;
const dimensionedVector F("F", ((2.0*Omega) & gHat)*gHat);
