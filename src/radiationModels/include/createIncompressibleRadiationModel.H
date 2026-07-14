    autoPtr<radiationModel> radiation(radiationModel::New(T));

    dimensionedScalar rhoCpRef
    (
        "rhoCpRef",
        dimensions::density*dimensions::specificHeatCapacity
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
            dimensions::density,
            physicalProperties
        );

        dimensionedScalar CpRef
        (
            "CpRef",
            dimensions::specificHeatCapacity,
            physicalProperties
        );

        rhoCpRef = rhoRef*CpRef;
    }
