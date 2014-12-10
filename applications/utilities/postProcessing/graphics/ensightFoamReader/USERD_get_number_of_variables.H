//======================================================================
int USERD_get_number_of_variables(void)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_number_of_variables" << endl << flush;
#endif

    return Num_variables;
}
