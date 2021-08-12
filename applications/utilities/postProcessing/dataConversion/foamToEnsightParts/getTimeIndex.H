// Read time index from */uniform/time, but treat 0 and constant specially
// or simply increment from the '-index' option if it was supplied

    label timeIndex = 0;

    if (optIndex)
    {
        timeIndex = indexingNumber++;
    }
    else if
    (
        runTime.timeName() != runTime.constant()
     && runTime.timeName() != "0"
    )
    {
        typeIOobject<IOdictionary> io
        (
            "time",
            runTime.timeName(),
            "uniform",
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        if (io.headerOk())
        {
            io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
            IOdictionary timeObject(io);

            timeObject.lookup("index") >> timeIndex;
        }
        else
        {
            Info<< "skip ... missing entry "
                << io.objectPath() << endl;
            continue;
        }
    }

    timeIndices.insert(timeIndex, timeDirs[timeI].value());
    Info<< "\nTime [" << timeIndex << "] = " << runTime.timeName() << nl;
