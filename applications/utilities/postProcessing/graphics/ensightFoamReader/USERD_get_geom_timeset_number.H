int USERD_get_geom_timeset_number
(
    void
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_geom_timeset_number" << endl
        << flush;
#endif

    Geom_timeset_number = 1;

    return Geom_timeset_number;
}
