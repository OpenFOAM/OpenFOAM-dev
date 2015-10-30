    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(T)
    );

    dimensionedScalar rhoCpRef
    (
        "rhoCpRef",
        dimDensity*dimEnergy/dimMass/dimTemperature,
        1.0
    );

    if (radiation->radiation())
    {
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false // Do not register
            )
        );

        dimensionedScalar rhoRef
        (
            "rhoRef",
            dimDensity,
            transportProperties
        );

        dimensionedScalar CpRef
        (
            "CpRef",
            dimSpecificHeatCapacity,
            transportProperties
        );

        rhoCpRef = rhoRef*CpRef;
    }
