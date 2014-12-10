    #include "readGravitationalAcceleration.H"

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
        mesh,
        dimensionedVector("U", dimVelocity, vector::zero)
    );

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phi", dimArea*dimVelocity, 0)
    );

    multiphaseSystem fluid(U, phi);

    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();
        const volScalarField& alpha = phase;

        U += alpha*phase.U();
        phi += fvc::interpolate(alpha)*phase.phi();
    }

    scalar slamDampCoeff
    (
        fluid.lookupOrDefault<scalar>("slamDampCoeff", 1)
    );

    dimensionedScalar maxSlamVelocity
    (
        "maxSlamVelocity",
        dimVelocity,
        fluid.lookupOrDefault<scalar>("maxSlamVelocity", GREAT)
    );

    // dimensionedScalar pMin
    // (
    //     "pMin",
    //     dimPressure,
    //     fluid.lookup("pMin")
    // );

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
        fluid.rho()
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


    autoPtr<incompressible::LESModel> sgsModel
    (
        incompressible::LESModel::New(U, phi, fluid)
    );
