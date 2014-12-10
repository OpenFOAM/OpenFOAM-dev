//======================================================================
int USERD_get_descrip_lines
(
    int which_type,
    int which_var,
    int imag_data,
    char line1[Z_BUFL],
    char line2[Z_BUFL]
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_descrip_lines" << endl
        << flush;
#endif

    if (which_type == Z_GEOM)
    {
        strncpy(line1, meshName, Z_BUFL);
        strncpy(line2, "", Z_BUFL);
    }
    else
    {
        strncpy(line1, "WHERE IS THIS LINE USED I WONDER???", Z_BUFL);
    }

#ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_descrip_lines" << endl
        << flush;
#endif
    return Z_OK;
}

