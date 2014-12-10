//======================================================================
//  Setting name in the gui, and specifying one or two input fields
//======================================================================
int USERD_get_name_of_reader
(
    char reader_name[Z_MAX_USERD_NAME],
    int *two_fields
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_name_of_reader" << endl << flush;
#endif

    strncpy(reader_name, readerName, Z_MAX_USERD_NAME);
    *two_fields = FALSE;

    return Z_OK;
}
