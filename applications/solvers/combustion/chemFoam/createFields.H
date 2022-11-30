    if (mesh.nCells() != 1)
    {
        FatalErrorInFunction
            << "Solver only applicable to single cell cases"
            << exit(FatalError);
    }

    Info<< "Reading initial conditions.\n" << endl;
    IOdictionary initialConditions
    (
        IOobject
        (
            "initialConditions",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar p0 = initialConditions.lookup<scalar>("p");
    scalar T0 = initialConditions.lookup<scalar>("T");

    #include "createBaseFields.H"

    Info<< "Reading thermophysical properties\n" << endl;
    autoPtr<fluidMulticomponentThermo> pThermo
    (
        fluidMulticomponentThermo::New(mesh)
    );
    fluidMulticomponentThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h");

    autoPtr<basicChemistryModel> pChemistry(basicChemistryModel::New(thermo));

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.renameRho()
    );

    volScalarField Rspecific
    (
        IOobject
        (
            "Rspecific",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimensionSet(dimEnergy/dimMass/dimTemperature),
            0
        )
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero)
    );

    #include "createPhi.H"

    OFstream post(args.path()/"chemFoam.out");
    post<< "# Time" << token::TAB << "Temperature [K]" << token::TAB
        << "Pressure [Pa]" << endl;
