
int USERD_get_number_of_model_parts(void)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_number_of_model_parts" << endl << flush;
#endif

    return Numparts_available;
}
