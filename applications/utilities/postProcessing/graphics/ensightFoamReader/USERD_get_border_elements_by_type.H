
// Not called if USERD_border_availability returns Z_ERR

int USERD_get_border_elements_by_type
(
    int part_number,
    int element_type,
    int **conn_array,
    short *parent_element_type,
    int *parnet_element_type
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_border_elements_by_type" << endl
        << flush;
#endif

    return Z_ERR;
}
