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
                false  // do not register!
            )
        );

        dimensionedScalar rhoRef(transportProperties.lookup("rhoRef"));
        dimensionedScalar CpRef(transportProperties.lookup("CpRef"));

        rhoCpRef = rhoRef*CpRef;
    }
