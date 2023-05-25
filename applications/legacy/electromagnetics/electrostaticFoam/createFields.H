    Info<< "Reading physicalProperties\n" << endl;

    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar epsilon0
    (
        physicalProperties.lookup("epsilon0")
    );

    dimensionedScalar k
    (
        physicalProperties.lookup("k")
    );


    Info<< "Reading field phi\n" << endl;
    volScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Calculating field rhoFlux\n" << endl;
    surfaceScalarField rhoFlux
    (
        IOobject
        (
            "rhoFlux",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -k*mesh.magSf()*fvc::snGrad(phi)
    );
