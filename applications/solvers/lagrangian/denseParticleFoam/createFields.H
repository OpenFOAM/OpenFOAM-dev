#include "readGravitationalAcceleration.H"

word continuousPhaseName
(
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    ).lookup("continuousPhaseName")
);

Info<< "Reading field U\n" << endl;
volVectorField Uc
(
    IOobject
    (
        IOobject::groupName("U", continuousPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading/calculating continuous-phase face flux field phic\n"
    << endl;

surfaceScalarField phic
(
    IOobject
    (
        IOobject::groupName("phi", continuousPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Uc) & mesh.Sf()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;

singlePhaseTransportModel continuousPhaseTransport(Uc, phic);

dimensionedScalar rhocValue
(
    IOobject::groupName("rho", continuousPhaseName),
    dimDensity,
    continuousPhaseTransport.lookup
    (
        IOobject::groupName("rho", continuousPhaseName)
    )
);

volScalarField rhoc
(
    IOobject
    (
        rhocValue.name(),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhocValue
);

volScalarField muc
(
    IOobject
    (
        IOobject::groupName("mu", continuousPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhoc*continuousPhaseTransport.nu()
);

Info << "Creating field alphac\n" << endl;
// alphac must be constructed before the cloud
// so that the drag-models can find it
volScalarField alphac
(
    IOobject
    (
        IOobject::groupName("alpha", continuousPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0)
);

Info<< "Constructing clouds" << endl;
parcelCloudList clouds(rhoc, Uc, muc, g);

// Particle fraction upper limit
scalar alphacMin
(
    1 - mesh.solverDict(alphac.name()).lookup<scalar>("max")
);

// Update alphac from the particle locations
alphac = max(1.0 - clouds.theta(), alphacMin);
alphac.correctBoundaryConditions();

surfaceScalarField alphacf("alphacf", fvc::interpolate(alphac));

surfaceScalarField alphaPhic
(
    IOobject::groupName
    (
        "alphaPhi",
        continuousPhaseName
    ),
    alphacf*phic
);

autoPtr<PhaseIncompressibleMomentumTransportModel<singlePhaseTransportModel>>
continuousPhaseTurbulence
(
    PhaseIncompressibleMomentumTransportModel<singlePhaseTransportModel>::New
    (
        alphac,
        Uc,
        alphaPhic,
        phic,
        continuousPhaseTransport
    )
);

#include "createFvOptions.H"
