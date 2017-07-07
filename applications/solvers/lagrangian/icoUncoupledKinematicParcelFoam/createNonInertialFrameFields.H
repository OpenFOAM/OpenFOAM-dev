    Info<< "Reading non-inertial frame fields" << endl;

    IOobject linearAccelerationHeader
    (
        "linearAcceleration",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> linearAccelerationPtr;

    if
    (
        linearAccelerationHeader.typeHeaderOk<uniformDimensionedVectorField>
        (
            true
        )
    )
    {
        Info<< "    Reading " << linearAccelerationHeader.name() << endl;

        linearAccelerationPtr.reset
        (
            new uniformDimensionedVectorField(linearAccelerationHeader)
        );
    }


    IOobject angularVelocityHeader
    (
        "angularVelocity",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> angularVelocityPtr;

    if (angularVelocityHeader.typeHeaderOk<uniformDimensionedVectorField>(true))
    {
        Info<< "    Reading " << angularVelocityHeader.name() << endl;

        angularVelocityPtr.reset
        (
            new uniformDimensionedVectorField(angularVelocityHeader)
        );
    }


    IOobject angularAccelerationHeader
    (
        "angularAcceleration",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> angularAccelerationPtr;

    if
    (
        angularAccelerationHeader.typeHeaderOk<uniformDimensionedVectorField>
        (
            true
        )
    )
    {
        Info<< "    Reading " << angularAccelerationHeader.name() << endl;

        angularAccelerationPtr.reset
        (
            new uniformDimensionedVectorField(angularAccelerationHeader)
        );
    }


    IOobject centreOfRotationHeader
    (
        "centreOfRotation",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> centreOfRotationPtr;

    if
    (
        centreOfRotationHeader.typeHeaderOk<uniformDimensionedVectorField>
        (
            true
        )
    )
    {
        Info<< "    Reading " << centreOfRotationHeader.name() << endl;

        centreOfRotationPtr.reset
        (
            new uniformDimensionedVectorField(centreOfRotationHeader)
        );
    }
