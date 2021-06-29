Info<< "Create engine time\n" << endl;

autoPtr<engineTime> runTimePtr
(
    engineTime::New
    (
        Time::controlDictName,
        args.rootPath(),
        args.globalCaseName(),
        args.caseName()
    )
);

engineTime& runTime = runTimePtr();
