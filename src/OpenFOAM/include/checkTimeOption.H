// Check -time and -latestTime options

if (args.optionFound("time"))
{
    Foam::scalar timeValue = args.optionRead<scalar>("time");

    startTime = Foam::Time::findClosestTimeIndex(Times, timeValue);
}

if (args.optionFound("latestTime"))
{
    startTime = Times.size() - 1;
}
