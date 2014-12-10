    Info<< "Reading thermodynamicProperties\n" << endl;

    IOdictionary thermodynamicProperties
    (
        IOobject
        (
            "thermodynamicProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rho0
    (
        thermodynamicProperties.lookup("rho0")
    );

    dimensionedScalar p0
    (
        thermodynamicProperties.lookup("p0")
    );

    dimensionedScalar psi
    (
        thermodynamicProperties.lookup("psi")
    );

    // Density offset, i.e. the constant part of the density
    dimensionedScalar rhoO("rhoO", rho0 - psi*p0);
