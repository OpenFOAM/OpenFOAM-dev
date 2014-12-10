int USERD_get_number_of_timesets
(
    void
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_number_of_timesets" << endl
        << flush;
#endif

    Num_timesets = 1;
    return  Num_timesets;
}
