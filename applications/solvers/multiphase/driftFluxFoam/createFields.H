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


    Info<< "Reading transportProperties\n" << endl;
    incompressibleTwoPhaseInteractingMixture mixture(U, phi);

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    const dimensionedScalar& rho1 = mixture.rhod();
    const dimensionedScalar& rho2 = mixture.rhoc();

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Mixture density
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture.rho()
    );

    // Mass flux
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );


    // Relative Velocity
    autoPtr<relativeVelocityModel> UdmModelPtr
    (
        relativeVelocityModel::New
        (
            transportProperties,
            mixture
        )
    );

    relativeVelocityModel& UdmModel(UdmModelPtr());


    // Construct compressible turbulence model
    autoPtr
    <
        CompressibleTurbulenceModel<incompressibleTwoPhaseInteractingMixture>
    > turbulence
    (
        CompressibleTurbulenceModel<incompressibleTwoPhaseInteractingMixture>
        ::New(rho, U, rhoPhi, mixture)
    );

    Info<< "Calculating field (g.h)f\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("gh", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }


    // MULES Correction
    tmp<surfaceScalarField> tphiAlphaCorr0;
