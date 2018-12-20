    Info<< "\nReading pRef" << endl;
    uniformDimensionedScalarField pRef
    (
        IOobject
        (
            "pRef",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimPressure, 0)
    );
