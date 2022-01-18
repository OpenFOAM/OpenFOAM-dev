{
    instantList sourceTimes = runTimeSource.times();

    if (sourceTimes.empty())
    {
        FatalErrorInFunction << "No result times in source "
            << runTimeSource.caseName()
            << exit(FatalError);
    }

    label sourceTimeIndex = runTimeSource.timeIndex();
    if (args.optionFound("sourceTime"))
    {
        if (args["sourceTime"] == "latestTime")
        {
            sourceTimeIndex = sourceTimes.size() - 1;
        }
        else
        {
            sourceTimeIndex = Time::findClosestTimeIndex
            (
                sourceTimes,
                args.optionRead<scalar>("sourceTime")
            );
        }
    }
    else
    {
        sourceTimeIndex = Time::findClosestTimeIndex
        (
            sourceTimes,
            runTimeTarget.time().userTimeValue()
        );
    }

    runTimeSource.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);

    Info<< "\nSource time: " << runTimeSource.userTimeName()
        << "\nTarget time: " << runTimeTarget.userTimeName()
        << endl;
}
