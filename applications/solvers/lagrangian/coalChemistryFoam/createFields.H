    Info<< "Creating combustion model\n" << endl;

    autoPtr<combustionModels::psiCombustionModel> combustion
    (
        combustionModels::psiCombustionModel::New(mesh)
    );

    psiReactionThermo& thermo = combustion->thermo();
    thermo.validate(args.executable(), "h", "e");

    SLGThermo slgThermo(mesh, thermo);

    basicMultiComponentMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    const word inertSpecie(thermo.lookup("inertSpecie"));

    if (!composition.contains(inertSpecie))
    {
        FatalErrorIn(args.executable())
            << "Specified inert specie '" << inertSpecie << "' not found in "
            << "species list. Available species:" << composition.species()
            << exit(FatalError);
    }

    volScalarField& p = thermo.p();
    const volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(thermo.he());

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    // lagrangian effective density field - used externally (optional)
    volScalarField rhoEffLagrangian
    (
        IOobject
        (
            "rhoEffLagrangian",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimDensity, 0.0)
    );

    // dynamic pressure field - used externally (optional)
    volScalarField pDyn
    (
        IOobject
        (
            "pDyn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    );


    Info<< "\nReading field U\n" << endl;
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

    #include "compressibleCreatePhi.H"

    Info<< "Creating turbulence model\n" << endl;
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

    // Set the turbulence into the combustion model
    combustion->setTurbulence(turbulence());

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
    volScalarField K("K", 0.5*magSqr(U));

    volScalarField dQ
    (
        IOobject
        (
            "dQ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
    );
