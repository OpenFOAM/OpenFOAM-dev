//======================================================================
void USERD_set_time_set_and_step
(
    int timeset_number,
    int time_step
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_set_time_set_and_step" << endl << flush;
#endif
    // update the global pointers and variables
    // to the current time-step

    // at exit time_step < 0
    if (time_step >= 0)
    {
        Time& runTime = *runTimePtr;
        Current_time_step = time_step;
        // add 1, since the first timestep is 'constant'

        if (time_step == 0)
        {
            runTime.setTime
            (
                timeDirs[Current_time_step],
                Current_time_step
            );
        }
        else
        {
            runTime.setTime
            (
                timeDirs[Current_time_step + 1],
                Current_time_step + 1
            );
        }

        meshPtr->readUpdate();

        if (time_step == 0)
        {
            runTime.setTime
            (
                timeDirs[Current_time_step + 1],
                Current_time_step + 1
            );
        }

        if (Numparts_available > nPatches+1)
        {
            delete sprayPtr;
            sprayPtr =
                new lagrangian::Cloud<passiveParticle>
                (
                    *meshPtr,
                    lagrangian::cloud::defaultName
                );
        }
    }
}
