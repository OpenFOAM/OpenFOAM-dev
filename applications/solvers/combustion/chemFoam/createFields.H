    if (mesh.nCells() != 1)
    {
        FatalErrorIn(args.executable())
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
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar p0 = readScalar(initialConditions.lookup("p"));
    scalar T0 = readScalar(initialConditions.lookup("T"));

    #include "createBaseFields.H"

    Info<< nl << "Reading thermophysicalProperties" << endl;
    autoPtr<psiChemistryModel> pChemistry(psiChemistryModel::New(mesh));

    psiChemistryModel& chemistry = pChemistry();
    scalar dtChem = refCast<const psiChemistryModel>(chemistry).deltaTChem()[0];

    psiReactionThermo& thermo = chemistry.thermo();
    thermo.validate(args.executable(), "h");

    basicMultiComponentMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    volScalarField& p = thermo.p();

    volScalarField Rspecific
    (
        IOobject
        (
            "Rspecific",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(dimEnergy/dimMass/dimTemperature),
            0.0
        )
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero),
        p.boundaryField().types()
    );

    #include "createPhi.H"

    Info << "Creating turbulence model.\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    OFstream post(args.path()/"chemFoam.out");
    post<< "# Time" << token::TAB << "Temperature [K]" << token::TAB
        << "Pressure [Pa]" << endl;


