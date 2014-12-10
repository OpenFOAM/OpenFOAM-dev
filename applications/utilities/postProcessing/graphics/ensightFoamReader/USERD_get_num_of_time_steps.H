//======================================================================
// don't use multiple time sets...NN
//======================================================================
int USERD_get_num_of_time_steps
(
    int timeset_number
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_num_of_time_steps" << endl
        << flush;
#endif

    return Num_time_steps;
}
