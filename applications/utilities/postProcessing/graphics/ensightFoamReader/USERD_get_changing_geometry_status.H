//======================================================================
int USERD_get_changing_geometry_status(void)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_changing_geometry_status" << endl << flush;
#endif

    // Choose the most general option
    return Z_CHANGE_CONN;
}

