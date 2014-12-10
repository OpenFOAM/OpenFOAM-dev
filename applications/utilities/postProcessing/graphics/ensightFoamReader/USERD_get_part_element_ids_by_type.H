int USERD_get_part_element_ids_by_type
(
    int part_number,
    int element_type,
    int *elemid_array
)
{

#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_part_element_ids_by_type" << endl
        << "part_number = " << part_number << endl
        << "element_type = " << element_type << endl << flush;
#endif

    if (part_number == 1)
    {
        const cellShapeList& cellShapes = meshPtr->cellShapes();
        const cellList& cells = meshPtr->cells();

        label nCells = cells.size();

        label nPen06 = 0;
        label nHex08 = 0;
        label nPyr05 = 0;
        label nTet04 = 0;
        label nFaced = 0;

        if (element_type == Z_HEX08)
        {
            for (label n=0; n<nCells; n++)
            {
                label nFaces = cells[n].size();
                labelList points = cellShapes[n];

                if ((nFaces == 6) && (points.size() == 8))
                {
                    elemid_array[nHex08++] = n + 1;
                }
            }
        }
        else if (element_type == Z_PEN06)
        {
            for (label n=0; n<nCells; n++)
            {
                label nFaces = cells[n].size();
                labelList points = cellShapes[n];

                if ((nFaces == 5) && (points.size() == 6))
                {
                    elemid_array[nPen06++] = n + 1;
                }
            }
        }
        else if (element_type == Z_PYR05)
        {
            for (label n=0; n<nCells; n++)
            {
                label nFaces = cells[n].size();
                labelList points = cellShapes[n];

                if ((nFaces == 5) && (points.size() == 5))
                {
                    elemid_array[nPyr05++] = n + 1;
                }
            }
        }
        else if (element_type == Z_TET04)
        {
            for (label n=0; n<nCells; n++)
            {
                label nFaces = cells[n].size();
                labelList points = cellShapes[n];

                if ((nFaces == 4) && (points.size() == 4))
                {
                    elemid_array[nTet04++] = n + 1;
                }
            }
        }
        else if (element_type == Z_NFACED)
        {
            for (label n=0; n<nCells; n++)
            {
                label nFaces = cells[n].size();
                labelList points = cellShapes[n];
                if ((nFaces == 6) && (points.size() == 8))
                {}
                else if ((nFaces == 5) && (points.size() == 6))
                {}
                else if ((nFaces == 5) && (points.size() == 5))
                {}
                else if ((nFaces == 4) && (points.size() == 4))
                {}
                else
                {
                    elemid_array[nFaced++] = n + 1;
                }

            }
        }
    }
    else if (part_number < nPatches+2)
    {

        const polyBoundaryMesh& bMesh = meshPtr->boundaryMesh();
        label patchi = part_number - 2;

        label nTri03 = 0;
        label nQuad04 = 0;
        label nPoly = 0;

        if (element_type == Z_TRI03)
        {
            forAll(bMesh[patchi], facei)
            {
                if (bMesh[patchi][facei].size() == 3)
                {
                    elemid_array[nTri03++] = facei + 1;
                }
            }
        }
        else if (element_type == Z_QUA04)
        {
            forAll(bMesh[patchi], facei)
            {
                if (bMesh[patchi][facei].size() == 4)
                {
                    elemid_array[nQuad04++] = facei + 1;
                }
            }
        }
        else if (element_type == Z_NSIDED)
        {
            forAll(bMesh[patchi], facei)
            {
                label nPoints = bMesh[patchi][facei].size();
                if ((nPoints != 3) && (nPoints != 4))
                {
                    elemid_array[nPoly++] = facei + 1;
                }
            }
        }

    }
    else if (part_number == nPatches+2)
    {
        forAll(*sprayPtr, n)
        {
            elemid_array[n] = n + 1;
        }
    }
    else
    {
        return Z_ERR;
    }

#ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_part_element_ids_by_type" << endl << flush;
#endif

    return Z_OK;
}
