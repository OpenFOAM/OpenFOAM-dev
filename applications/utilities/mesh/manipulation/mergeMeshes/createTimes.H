    Info<< nl << "Create Times" << endl;

    const fileName masterCasePath = masterCase.path();
    const fileName masterCaseName = masterCase.name();

    Time runTimeMaster
    (
        Time::controlDictName,
        masterCasePath,
        masterCaseName
    );

    const fileName addCasePath = addCase.path();
    const fileName addCaseName = addCase.name();

    Time runTimeToAdd
    (
        Time::controlDictName,
        addCasePath,
        addCaseName
    );

