//
// setConstantMeshDictionaryIO.H
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    fileName dictPath = "";
    if (args.optionFound("dict"))
    {
        dictPath = args["dict"];
        if (isDir(dictPath))
        {
            dictPath = dictPath / dictName;
        }
    }

    IOobject dictIO
    (
        dictName,
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );

    if (dictPath.size())
    {
        dictIO = IOobject
        (
            dictPath,
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        );
    }
