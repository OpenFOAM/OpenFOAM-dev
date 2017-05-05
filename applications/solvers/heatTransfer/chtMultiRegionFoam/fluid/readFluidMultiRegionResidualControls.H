    const dictionary& residualControl =
        mesh.solutionDict().subDict("PIMPLE").subOrEmptyDict("residualControl");

    scalar UTol = -1.;
    scalar ETol = -1.;
    scalar p_rghTol = -1.;

    if (!residualControl.empty())
    {
        if (!residualControl.subOrEmptyDict("U").empty())
        {
            UTol = readScalar(residualControl.subDict("U").lookup("tolerance"));
        }

        if (!residualControl.subOrEmptyDict("p_rgh").empty())
        {
            p_rghTol =
                readScalar
                (
                    residualControl.subDict("p_rgh").lookup("tolerance")
                );
        }

        if (!residualControl.subOrEmptyDict("h").empty())
        {
            ETol = readScalar(residualControl.subDict("h").lookup("tolerance"));
        }

        // Residual control used?
        if (UTol != -1 || ETol != -1 || p_rghTol != -1)
        {
            residualControlUsed = true;
            resReachedFluid = false;
        }
    }
