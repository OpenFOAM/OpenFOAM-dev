//======================================================================
// this is based on the current time step.
//======================================================================
int USERD_get_gold_part_build_info
(
    int *part_numbers,
    int *part_types,
    char *part_descriptions[Z_BUFL],
    int *number_of_nodes,
    int *number_of_elements[Z_MAXTYPE],
    int *ijk_dimensions[9],
    int *iblanking_options[6]
)
{
    #ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_gold_part_build_info" << endl << flush;
    #endif

    const cellShapeList& cellShapes = meshPtr->cellShapes();
    const cellList& cells = meshPtr->cells();

    label nCells = cells.size();

    // all parts are unstructured
    for (label n = 0; n<Numparts_available; n++)
    {
        part_numbers[n] = n + 1;
        part_types[n]   = Z_UNSTRUCTURED;
    }

    strncpy(part_descriptions[0], meshName, Z_BUFL);

    for (label i=0; i<nPatches; i++)
    {
        word patchName(meshPtr->boundary()[i].name());
        strncpy(part_descriptions[i+1], patchName.c_str(), Z_BUFL);
    }

    label nHex08 = 0;
    label nPen06 = 0;
    label nPyr05 = 0;
    label nTet04 = 0;
    label nFaced = 0;

    for (label n=0; n<nCells; n++)
    {
        label nFacesInCell = cells[n].size();
        labelList points = cellShapes[n];

        if ((nFacesInCell == 6) && (points.size() == 8))
        {
            nHex08++;
        }
        else if ((nFacesInCell == 4) && (points.size() == 4))
        {
            nTet04++;
        }
        else if (nFacesInCell == 5)
        {
            if (points.size() == 6)
            {
                nPen06++;
            }
            else if (points.size() == 5)
            {
                nPyr05++;
            }
            else
            {
                nFaced++;
            }
        }
        else
        {
            nFaced++;
        }
    }

    for (label n=0; n < Z_MAXTYPE; n++)
    {
        for (label i=0; i<Numparts_available; i++)
        {
            number_of_elements[i][n] = 0;
        }
    }

    number_of_elements[0][Z_TET04] = nTet04;
    number_of_elements[0][Z_PYR05] = nPyr05;
    number_of_elements[0][Z_HEX08] = nHex08;
    number_of_elements[0][Z_PEN06] = nPen06;
    number_of_elements[0][Z_NFACED] = nFaced;

    /*
    Info<< "nTet04 = " << nTet04 << endl;
    Info<< "nPyr05 = " << nPyr05 << endl;
    Info<< "nHex08 = " << nHex08 << endl;
    Info<< "nPen06 = " << nPen06 << endl;
    Info<< "nFaced = " << nFaced << endl;
        */

    number_of_nodes[0] = meshPtr->nPoints();

    const polyBoundaryMesh& bMesh = meshPtr->boundaryMesh();

    for (label i=0; i<nPatches; i++)
    {
        label nTri03 = 0;
        label nQuad04 = 0;
        label nPoly = 0;

        forAll(bMesh[i], n)
        {
            label nPoints = bMesh[i][n].size();

            if (nPoints == 3)
            {
                nTri03++;
            }
            else  if (nPoints == 4)
            {
                nQuad04++;
            }
            else
            {
                nPoly++;
            }
        }

        number_of_elements[i+1][Z_TRI03] = nTri03;
        number_of_elements[i+1][Z_QUA04] = nQuad04;
        number_of_elements[i+1][Z_NSIDED] = nPoly;

        number_of_nodes[i+1] = bMesh[i].points().size();
    }

    if (Numparts_available > nPatches+1)
    {
        strncpy
        (
            part_descriptions[nPatches+1],
            lagrangian::cloud::prefix.c_str(),
            Z_BUFL
        );
        number_of_elements[nPatches+1][Z_POINT] = sprayPtr->size();
        number_of_nodes[nPatches+1] = sprayPtr->size();
    }

#ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_gold_part_build_info" << endl << flush;
#endif

    return Z_OK;
}
