    if (runTime.controlDict().lookupOrDefault("suppressSolverInfo", false))
    {
        lduMatrix::debug = 0;
    }

    Switch adjustTimeStep(runTime.controlDict().lookup("adjustTimeStep"));

    scalar maxDeltaT(readScalar(runTime.controlDict().lookup("maxDeltaT")));
