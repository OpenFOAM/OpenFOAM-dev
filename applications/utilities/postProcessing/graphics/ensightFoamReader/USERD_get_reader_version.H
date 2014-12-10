int USERD_get_reader_version
(
    char version_number[Z_MAX_USERD_NAME]
)
{

#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_reader_version" << endl;
#endif

    strncpy(version_number, readerVersion, Z_MAX_USERD_NAME);

#ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_reader_version" << endl;
#endif

    return Z_OK;
}
