// Read time index from */uniform/time, but treat 0 and constant specially

    word timeName = "0";

    if
    (
        runTime.timeName() != runTime.constant()
     && runTime.timeName() != "0"
    )
    {
        IOobject io
        (
            "time",
            runTime.timeName(),
            "uniform",
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        if (io.typeHeaderOk<IOdictionary>(true))
        {
            IOdictionary timeObject
            (
                IOobject
                (
                    "time",
                    runTime.timeName(),
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
            timeName = runTime.timeName();
        }
    }

    Info<< "\nTime [" << timeName << "] = " << runTime.timeName() << nl;
