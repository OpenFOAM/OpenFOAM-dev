    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<psiThermo> pThermo
    (
        psiThermo::New(mesh)
    );
    psiThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
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

    #include "compressibleCreatePhi.H"

    dimensionedScalar rhoMax(pimple.dict().lookup("rhoMax"));
    dimensionedScalar rhoMin(pimple.dict().lookup("rhoMin"));

    Info<< "Creating turbulence model\n" << endl;

    autoPtr<CompressibleTurbulenceModel<fluidThermo> >
    turbulence
    (
        CompressibleTurbulenceModel<fluidThermo>::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

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
