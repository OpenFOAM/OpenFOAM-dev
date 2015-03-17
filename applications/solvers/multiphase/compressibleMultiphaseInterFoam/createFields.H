    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
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

    #include "createPhi.H"

    Info<< "Constructing multiphaseMixtureThermo\n" << endl;
    multiphaseMixtureThermo mixture(U, phi);

    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        mixture.rho()
    );

    dimensionedScalar pMin(mixture.lookup("pMin"));


    #include "readGravitationalAcceleration.H"
    #include "readhRef.H"

    Info<< "Calculating field g.h\n" << endl;
    dimensionedScalar ghRef(g & (cmptMag(g.value())/mag(g.value()))*hRef);
    volScalarField gh("gh", (g & mesh.C()) - ghRef);
    surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);


    // Construct compressible turbulence model
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            mixture.rhoPhi(),
            mixture
        )
    );

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));
