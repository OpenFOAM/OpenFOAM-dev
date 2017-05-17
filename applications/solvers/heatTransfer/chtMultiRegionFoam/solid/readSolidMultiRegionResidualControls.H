    const dictionary& residualControl =
        mesh.solutionDict().subDict("PIMPLE").subOrEmptyDict("residualControl");

    scalar hTol = -1.;

    if (!residualControl.empty())
    {
        if (!residualControl.subOrEmptyDict("h").empty())
        {
            hTol = readScalar(residualControl.subDict("h").lookup("tolerance"));

            // Used residual control for actual solid region
            if (hTol != -1)
            {
                residualControlUsed = true;
                resReachedSolid = false;
            }
        }
    }
