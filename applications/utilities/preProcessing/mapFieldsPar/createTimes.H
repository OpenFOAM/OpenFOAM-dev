    Info<< "\nCreate databases as time" << endl;

    HashTable<string> srcOptions(args.options());
    srcOptions.erase("case");
    srcOptions.insert("case", fileName(rootDirSource/caseDirSource));

    argList argsSrc(args, srcOptions, false, false, false);
    if (!argsSrc.checkRootCase())
    {
        FatalError.exit();
    }

    Time runTimeSource(Time::controlDictName, argsSrc);

    Time runTimeTarget(Time::controlDictName, args);
