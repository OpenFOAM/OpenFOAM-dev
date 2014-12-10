    word dictName("noiseDict");
    if (args.optionFound("dict"))
    {
        dictName = args["dict"];
    }

    IOdictionary dict
    (
        IOobject
        (
            dictName,
            runTime.system(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // reference pressure
    scalar pRef = dict.lookupOrDefault("pRef", 0.0);

    // number of samples in sampling window
    label N = dict.lookupOrDefault("N", 65536);

    // number of sampling windows
    label nw = dict.lookupOrDefault("nw", 100);

    // lower frequency of frequency band
    scalar f1 = dict.lookupOrDefault("f1", 25.0);

    // upper frequency of frequency band
    scalar fU = dict.lookupOrDefault("fU", 10000.0);

    // graph format
    word graphFormat = dict.lookupOrDefault<word>("graphFormat", "raw");
