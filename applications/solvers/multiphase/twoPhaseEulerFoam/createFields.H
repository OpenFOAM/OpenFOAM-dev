    Info<< "Creating twoPhaseSystem\n" << endl;

    twoPhaseSystem fluid(mesh, g);

    phaseModel& phase1 = fluid.phase1();
    phaseModel& phase2 = fluid.phase2();

    volScalarField& alpha1 = phase1;
    volScalarField& alpha2 = phase2;

    volVectorField& U1 = phase1.U();
    surfaceScalarField& phi1 = phase1.phi();
    surfaceScalarField& alphaPhi1 = phase1.alphaPhi();
    surfaceScalarField& alphaRhoPhi1 = phase1.alphaRhoPhi();

    volVectorField& U2 = phase2.U();
    surfaceScalarField& phi2 = phase2.phi();
    surfaceScalarField& alphaPhi2 = phase2.alphaPhi();
    surfaceScalarField& alphaRhoPhi2 = phase2.alphaRhoPhi();

    surfaceScalarField& phi = fluid.phi();

    dimensionedScalar pMin
    (
        "pMin",
        dimPressure,
        fluid.lookup("pMin")
    );

    rhoThermo& thermo1 = phase1.thermo();
    rhoThermo& thermo2 = phase2.thermo();

    volScalarField& p = thermo1.p();

    volScalarField& rho1 = thermo1.rho();
    const volScalarField& psi1 = thermo1.psi();

    volScalarField& rho2 = thermo2.rho();
    const volScalarField& psi2 = thermo2.psi();

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.U()
    );

    Info<< "Calculating field DDtU1 and DDtU2\n" << endl;

    volVectorField DDtU1
    (
        "DDtU1",
        fvc::ddt(U1)
      + fvc::div(phi1, U1)
      - fvc::div(phi1)*U1
    );

    volVectorField DDtU2
    (
        "DDtU2",
        fvc::ddt(U2)
      + fvc::div(phi2, U2)
      - fvc::div(phi2)*U2
    );

    volScalarField rAU1
    (
        IOobject
        (
            IOobject::groupName("rAU", phase1.name()),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, 3, 1, 0, 0), 0.0)
    );

    volScalarField rAU2
    (
        IOobject
        (
            IOobject::groupName("rAU", phase2.name()),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, 3, 1, 0, 0), 0.0)
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
    );


    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K1(IOobject::groupName("K", phase1.name()), 0.5*magSqr(U1));
    volScalarField K2(IOobject::groupName("K", phase2.name()), 0.5*magSqr(U2));
