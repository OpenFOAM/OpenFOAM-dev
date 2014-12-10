
int USERD_get_nfaced_nodes_per_face
(
    int part_number,
    int *nfaced_npf_array
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_nfaced_nodes_per_face"
        << ", part_number = " << part_number
        << endl
        << flush;
#endif

    if (part_number == 1)
    {
        const cellShapeList& cellShapes = meshPtr->cellShapes();
        const cellList& cells = meshPtr->cells();
        const faceList& faces = meshPtr->faces();

        label nCells = cellShapes.size();
        label nFaced = 0;
        for (label n=0; n<nCells; n++)
        {
            label nFacesInCell = cells[n].size();
            labelList points = cellShapes[n];
            label nPoints = points.size();

            if ((nFacesInCell == 6) && (nPoints == 8))
            {}
            else if ((nFacesInCell == 4) && (nPoints == 4))
            {}
            else if (nFacesInCell == 5)
            {
                if (nPoints == 6)
                {}
                else if (nPoints == 5)
                {}
                else
                {
                    for (label i=0; i<nFacesInCell; i++)
                    {
                        label facei = cells[n][i];
                        label nFacePoints = faces[facei].size();
                        nfaced_npf_array[nFaced++] = nFacePoints;
                    }
                }
            }
            else
            {
                for (label i=0; i<nFacesInCell; i++)
                {
                    label facei = cells[n][i];
                    label nFacePoints = faces[facei].size();
                    nfaced_npf_array[nFaced++] = nFacePoints;
                }
            }
        }

    }
    else if (part_number < nPatches+2)
    {
        return Z_ERR;
    }
    else
    {
        return Z_ERR;
    }
#ifdef ENSIGHTDEBUG
    Info<< "Exiting: USERD_get_nfaced_nodes_per_face" << endl
        << flush;
#endif
    return Z_OK;
}
