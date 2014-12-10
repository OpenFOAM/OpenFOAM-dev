//
// checkTimeOption.H
// ~~~~~~~~~~~~~~~~~
// check -time and -latestTime options

    if (args.optionFound("time"))
    {
        Foam::scalar timeValue = args.optionRead<scalar>("time");

        startTime = Foam::Time::findClosestTimeIndex(Times, timeValue);
        endTime = startTime + 1;
    }

    if (args.optionFound("latestTime"))
    {
        startTime = Times.size() - 1;
    }
