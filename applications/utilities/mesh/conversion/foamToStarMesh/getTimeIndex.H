// Read time index from */uniform/time, but treat 0 and constant specially

    word timeName = "0";

    if
    (
        runTime.name() != runTime.constant()
     && runTime.name() != "0"
    )
    {
        typeIOobject<IOdictionary> io
        (
            "time",
            runTime.name(),
            "uniform",
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        if (io.headerOk())
        {
            IOdictionary timeObject
            (
                IOobject
                (
                    "time",
                    runTime.name(),
                    "uniform",
                    runTime,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    false
                )
            );

            label index;
            timeObject.lookup("index") >> index;
            timeName = Foam::name(index);
        }
        else
        {
            timeName = runTime.name();
        }
    }

    Info<< "\nTime [" << timeName << "] = " << runTime.name() << nl;
