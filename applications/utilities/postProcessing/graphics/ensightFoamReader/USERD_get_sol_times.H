//======================================================================
// Negative values of the time is not allowed in EnSight.
// So for engines, where the time is CAD's we need to correct
// this so that all CAD's are positive. NN
//======================================================================
int USERD_get_sol_times
(
    int timeset_number,
    float *solution_times
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_sol_times\n" << timeDirs << endl;
#endif

    for (label n=0; n<Num_time_steps;n++)
    {
        solution_times[n] = timeDirs[n+1].value();
    }

    if (timeDirs[1].value() < 0)
    {
        scalar addCAD = 360.0;
        while (timeDirs[1].value() + addCAD < 0.0)
        {
            addCAD += 360.0;
        }
        for (label n=0; n<Num_time_steps;n++)
        {
            solution_times[n] += addCAD;

            Info<< "Time[" << n << "] = " << timeDirs[n+1].value()
                << " was corrected to " << solution_times[n]  << endl;
        }

    }

#ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_sol_times" << endl;
#endif

    return Z_OK;
}

