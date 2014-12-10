int USERD_get_nsided_conn
(
    int part_number,
    int *nsided_conn_array
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_nsided_conn"
        << ", part_number = " << part_number
        << endl
        << flush;
#endif
    if (part_number == 1)
    {
        Info<< "************* EEEEEEEEERRRRRRRRRRRRRRRRRR *************** "
            << endl << flush;

    }
    else if (part_number < nPatches+2)
    {
        //const cellList& cells = meshPtr->cells();
        //const faceList& faces = meshPtr->faces();

        label patchi = part_number - 2;
        const polyBoundaryMesh& bMesh = meshPtr->boundaryMesh();

        label np = 0;
        forAll(bMesh[patchi], facei)
        {
            label nPoints = bMesh[patchi][facei].size();
            if ((nPoints != 3) && (nPoints != 4))
            {
                for (label i=0; i<nPoints; i++)
                {
                    nsided_conn_array[np++] = bMesh[patchi][facei][i] + 1;
                }
            }
        }
    }
    else if (part_number == nPatches+2)
    {
        return Z_ERR;
    }
#ifdef ENSIGHTDEBUG
    Info<< "Exiting: USERD_get_nsided_conn" << endl
        << flush;
#endif
    return Z_OK;
}
