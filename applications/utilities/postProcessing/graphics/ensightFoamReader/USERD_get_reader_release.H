int USERD_get_reader_release
(
    char release_number[Z_MAX_USERD_NAME]
)
{

#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_reader_release" << endl;
#endif

    strncpy(release_number, Foam::FOAMbuild, Z_MAX_USERD_NAME);

#ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_reader_release" << endl;
#endif

    return Z_OK;
}
