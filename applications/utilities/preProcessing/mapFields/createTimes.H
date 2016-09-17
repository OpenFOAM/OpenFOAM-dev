    Info<< "\nCreate databases as time" << endl;

    const auto caseDirOrig = getEnv("FOAM_CASE");
    const auto caseNameOrig = getEnv("FOAM_CASE_NAME");
    setEnv("FOAM_CASE", rootDirSource/caseDirSource, true);
    setEnv("FOAM_CASE_NAME", caseDirSource, true);
    Time runTimeSource
    (
        Time::controlDictName,
        rootDirSource,
        caseDirSource
    );
    setEnv("FOAM_CASE", caseDirOrig, true);
    setEnv("FOAM_CASE_NAME", caseNameOrig, true);

    Time runTimeTarget
    (
        Time::controlDictName,
        rootDirTarget,
        caseDirTarget
    );
