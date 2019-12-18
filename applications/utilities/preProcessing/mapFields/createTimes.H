    Info<< "\nCreate databases as time" << endl;

    const string caseDirOrig = getEnv("FOAM_CASE");
    const string caseNameOrig = getEnv("FOAM_CASENAME");
    setEnv("FOAM_CASE", rootDirSource/caseDirSource, true);
    setEnv("FOAM_CASENAME", caseDirSource, true);
    Time runTimeSource
    (
        Time::controlDictName,
        rootDirSource,
        caseDirSource
    );
    setEnv("FOAM_CASE", caseDirOrig, true);
    setEnv("FOAM_CASENAME", caseNameOrig, true);

    Time runTimeTarget
    (
        Time::controlDictName,
        rootDirTarget,
        caseDirTarget
    );
