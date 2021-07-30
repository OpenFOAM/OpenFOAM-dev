    autoPtr<radiationModel> radiation(radiationModel::New(T));

    dimensionedScalar rhoCpRef
    (
        "rhoCpRef",
        dimDensity*dimEnergy/dimMass/dimTemperature,
        1.0
    );

    if (!isType<radiationModels::noRadiation>(radiation()))
    {
        IOdictionary physicalProperties
        (
            IOobject
            (
                "physicalProperties",
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
            physicalProperties
        );

        dimensionedScalar CpRef
        (
            "CpRef",
            dimSpecificHeatCapacity,
            physicalProperties
        );

        rhoCpRef = rhoRef*CpRef;
    }
