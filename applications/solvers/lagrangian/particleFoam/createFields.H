#include "readGravitationalAcceleration.H"

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));

dimensionedScalar rhoInfValue("rhoInf", dimDensity, viscosity);

volScalarField rhoInf
(
    IOobject
    (
        "rho",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhoInfValue
);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, viscosity)
);

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    viscosity->nu()*rhoInfValue
);

Info<< "Constructing clouds" << endl;
parcelClouds& clouds = parcelClouds::New(mesh, rhoInf, U, mu, g);

typeIOobject<volVectorField> Hheader
(
    "H",
    runTime.name(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
);

autoPtr<volVectorField> HPtr;

if (Hheader.headerOk())
{
    Info<< "\nReading field H\n" << endl;

    HPtr.reset(new volVectorField (Hheader, mesh));
}

typeIOobject<volVectorField> HdotGradHheader
(
    "HdotGradH",
    runTime.name(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
);

autoPtr<volVectorField> HdotGradHPtr;

if (HdotGradHheader.headerOk())
{
    Info<< "Reading field HdotGradH" << endl;

    HdotGradHPtr.reset(new volVectorField(HdotGradHheader, mesh));
}

#include "createNonInertialFrameFields.H"
