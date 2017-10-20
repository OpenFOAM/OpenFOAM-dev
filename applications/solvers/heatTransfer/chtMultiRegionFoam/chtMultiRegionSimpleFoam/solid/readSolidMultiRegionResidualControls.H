const dictionary& residualControl =
    mesh.solutionDict().subDict("SIMPLE").subOrEmptyDict("residualControl");

scalar hTol = -1;

if (!residualControl.empty())
{
    if (!residualControl.subOrEmptyDict("h").empty())
    {
        hTol = readScalar(residualControl.subDict("h").lookup("tolerance"));

        //- Used residual control for actual solid region
        if (hTol != -1)
        {
            residualControlUsed = true;
            resReachedSolid = false;
        }
    }

    if (firstIteration)
    {
        firstIteration = false;

        if (hTol != -1)
        {
            Info<< "SIMPLE convergence criteria\n";

            if (hTol != -1)
            {
                Info<< "   field h      tolerance " << hTol << "\n";
            }

            Info<< endl;
        }
    }
}
