//======================================================================
//
//======================================================================
int USERD_get_number_of_files_in_dataset(void)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_number_of_files_in_dataset" << endl << flush;
#endif

    // use 1 instead of 0 which gives an un-necessary warning.
    Num_dataset_files = 1;
    return Num_dataset_files;
}
