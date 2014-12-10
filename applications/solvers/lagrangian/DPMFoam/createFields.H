    Info<< "\nReading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    word contiuousPhaseName(transportProperties.lookup("contiuousPhaseName"));

    dimensionedScalar rhocValue
    (
        IOobject::groupName("rho", contiuousPhaseName),
        dimDensity,
        transportProperties.lookup
        (
            IOobject::groupName("rho", contiuousPhaseName)
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

    Info<< "Reading field U\n" << endl;
    volVectorField Uc
    (
        IOobject
        (
            IOobject::groupName("U", contiuousPhaseName),
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
            IOobject::groupName("phi", contiuousPhaseName),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(Uc) & mesh.Sf()
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

    Info<< "Creating turbulence model\n" << endl;

    singlePhaseTransportModel continuousPhaseTransport(Uc, phic);

    volScalarField muc
    (
        IOobject
        (
            IOobject::groupName("mu", contiuousPhaseName),
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
            IOobject::groupName("alpha", contiuousPhaseName),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    );

    word kinematicCloudName("kinematicCloud");
    args.optionReadIfPresent("cloudName", kinematicCloudName);

    Info<< "Constructing kinematicCloud " << kinematicCloudName << endl;
    basicKinematicTypeCloud kinematicCloud
    (
        kinematicCloudName,
        rhoc,
        Uc,
        muc,
        g
    );

    // Particle fraction upper limit
    scalar alphacMin
    (
        1.0
      - readScalar
        (
            kinematicCloud.particleProperties().subDict("constantProperties")
           .lookup("alphaMax")
        )
    );

    // Update alphac from the particle locations
    alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
    alphac.correctBoundaryConditions();

    surfaceScalarField alphacf("alphacf", fvc::interpolate(alphac));
    surfaceScalarField alphaPhic("alphaPhic", alphacf*phic);

    autoPtr<PhaseIncompressibleTurbulenceModel<singlePhaseTransportModel> >
    continuousPhaseTurbulence
    (
        PhaseIncompressibleTurbulenceModel<singlePhaseTransportModel>::New
        (
            alphac,
            Uc,
            alphaPhic,
            phic,
            continuousPhaseTransport
        )
    );
