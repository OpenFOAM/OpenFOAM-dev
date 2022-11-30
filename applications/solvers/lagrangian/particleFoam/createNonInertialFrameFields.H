    Info<< "Reading non-inertial frame fields" << endl;

    typeIOobject<uniformDimensionedVectorField> linearAccelerationHeader
    (
        "linearAcceleration",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> linearAccelerationPtr;

    if (linearAccelerationHeader.headerOk())
    {
        Info<< "    Reading " << linearAccelerationHeader.name() << endl;

        linearAccelerationPtr.reset
        (
            new uniformDimensionedVectorField(linearAccelerationHeader)
        );
    }


    typeIOobject<uniformDimensionedVectorField> angularVelocityHeader
    (
        "angularVelocity",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> angularVelocityPtr;

    if (angularVelocityHeader.headerOk())
    {
        Info<< "    Reading " << angularVelocityHeader.name() << endl;

        angularVelocityPtr.reset
        (
            new uniformDimensionedVectorField(angularVelocityHeader)
        );
    }


    typeIOobject<uniformDimensionedVectorField> angularAccelerationHeader
    (
        "angularAcceleration",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> angularAccelerationPtr;

    if (angularAccelerationHeader.headerOk())
    {
        Info<< "    Reading " << angularAccelerationHeader.name() << endl;

        angularAccelerationPtr.reset
        (
            new uniformDimensionedVectorField(angularAccelerationHeader)
        );
    }


    typeIOobject<uniformDimensionedVectorField> centreOfRotationHeader
    (
        "centreOfRotation",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> centreOfRotationPtr;

    if (centreOfRotationHeader.headerOk())
    {
        Info<< "    Reading " << centreOfRotationHeader.name() << endl;

        centreOfRotationPtr.reset
        (
            new uniformDimensionedVectorField(centreOfRotationHeader)
        );
    }
