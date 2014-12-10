
int USERD_get_nfaced_conn
(
    int part_number,
    int *nfaced_conn_array
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_nfaced_conn"
        << ", part_number = " << part_number
        << endl
        << flush;
#endif

    if (part_number == 1)
    {
        label nPoint = 0;
        const cellShapeList& cellShapes = meshPtr->cellShapes();
        const cellList& cells = meshPtr->cells();
        const faceList& faces = meshPtr->faces();
        label nCells = cellShapes.size();

        for (label n=0; n<nCells; n++)
        {
            label nFacesInCell = cells[n].size();
            labelList points = cellShapes[n];
            if ((nFacesInCell == 6) && (points.size() == 8))
            {}
            else if ((nFacesInCell == 4) && (points.size() == 4))
            {}
            else if (nFacesInCell == 5)
            {
                if (points.size() == 6)
                {}
                else if (points.size() == 5)
                {}
                else
                {
                    for (label i=0; i<nFacesInCell; i++)
                    {
                        label facei = cells[n][i];
                        label nPoints = faces[facei].size();
                        for (label j=0; j<nPoints; j++)
                        {
                            nfaced_conn_array[nPoint++] = faces[facei][j] + 1;
                        }
                    }
                }
            }
            else
            {
                for (label i=0; i<nFacesInCell; i++)
                {
                    label facei = cells[n][i];
                    label nPoints = faces[facei].size();
                    for (label j=0; j<nPoints; j++)
                    {
                        nfaced_conn_array[nPoint++] = faces[facei][j] + 1;
                    }
                }
            }
        }

    }
    else if (part_number < nPatches+2)
    {

    }
    else
    {
        return Z_ERR;
    }

#ifdef ENSIGHTDEBUG
    Info<< "Exiting: USERD_get_nfaced_conn" << endl
        << flush;
#endif
    return Z_OK;
}
