int USERD_get_part_elements_by_type
(
    int part_number,
    int element_type,
    int **conn_array
)
{
    #ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_part_elements_by_type" << nl
        << "part_number = " << part_number << nl
        << "element_type = " << element_type;
    if (element_type == Z_HEX08)
    {
        Info<< " Z_HEX08";
    }
    else if (element_type == Z_PEN06)
    {
        Info<< " Z_PEN06";
    }
    else if (element_type == Z_PYR05)
    {
        Info<< " Z_PYR05";
    }
    else if (element_type == Z_TET04)
    {
        Info<< " Z_TET04";
    }
    else if (element_type == Z_TRI03)
    {
        Info<< " Z_TRI03";
    }
    else if (element_type == Z_QUA04)
    {
        Info<< " Z_QUA04";
    }
    else if (element_type == Z_NFACED)
    {
        Info<< " Z_NFACED";
    }
    else if (element_type == Z_NSIDED)
    {
        Info<< " Z_NSIDED";
    }
    else
    {
        Info<< " unknown";
    }
    Info<< endl << flush;
    #endif

    if (part_number == 1)
    {
        const cellShapeList& cellShapes = meshPtr->cellShapes();

        //================================
        // hexahedron
        //================================
        if (element_type == Z_HEX08)
        {
            const cellModel& hex = *(cellModeller::lookup("hex"));

            label nHex08 = 0;
            forAll(cellShapes, celli)
            {
                const cellShape& cellShape = cellShapes[celli];
                const cellModel& cellModel = cellShape.model();

                if (cellModel == hex)
                {
                    forAll(cellShape, ip)
                    {
                        conn_array[nHex08][ip] = cellShape[ip] + 1;
                    }
                    nHex08++;
                }
            }
        }
        //================================
        // pentahedron
        //================================
        else if (element_type == Z_PEN06)
        {
            const cellModel& prism = *(cellModeller::lookup("prism"));

            label nPen06 = 0;
            forAll(cellShapes, celli)
            {
                const cellShape& cellShape = cellShapes[celli];
                const cellModel& cellModel = cellShape.model();

                if (cellModel == prism)
                {
                    forAll(cellShape, ip)
                    {
                        conn_array[nPen06][ip] = cellShape[ip] + 1;
                    }
                    nPen06++;
                }
            }
        }
        //================================
        // pyramid
        //================================
        else if (element_type == Z_PYR05)
        {
            const cellModel& pyr = *(cellModeller::lookup("pyr"));

            label nPyr05 = 0;
            forAll(cellShapes, celli)
            {
                const cellShape& cellShape = cellShapes[celli];
                const cellModel& cellModel = cellShape.model();

                if (cellModel == pyr)
                {
                    forAll(cellShape, ip)
                    {
                        conn_array[nPyr05][ip] = cellShape[ip] + 1;
                    }
                    nPyr05++;
                }
            }
        }
        //================================
        // tetrahedron
        //================================
        else if (element_type == Z_TET04)
        {
            const cellModel& tet = *(cellModeller::lookup("tet"));

            label nTet04 = 0;
            forAll(cellShapes, celli)
            {
                const cellShape& cellShape = cellShapes[celli];
                const cellModel& cellModel = cellShape.model();

                if (cellModel == tet)
                {
                    forAll(cellShape, ip)
                    {
                        conn_array[nTet04][ip] = cellShape[ip] + 1;
                    }
                    nTet04++;
                }
            }
        }
        //================================
        // polyhedra
        //================================
        else
        {
            label nCells = cellShapes.size();
            label nFaced = 0;
            const cellList cells = meshPtr->cells();

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
                        conn_array[nFaced++][0] = nFacesInCell;
                    }
                }
                else
                {
                    conn_array[nFaced++][0] = nFacesInCell;
                }
            }
        }
    }
    else if (part_number < nPatches+2)
    {
        label patchi = part_number - 2;
        const polyBoundaryMesh& bMesh = meshPtr->boundaryMesh();

        label nTri03 = 0;
        label nQuad04 = 0;
        if (element_type == Z_TRI03)
        {
            forAll(bMesh[patchi], n)
            {
                label nPoints = bMesh[patchi][n].size();
                if (nPoints == 3)
                {
                    for (label i=0; i<nPoints; i++)
                    {
                        label ip = bMesh[patchi][n][i];
                        conn_array[nTri03][i] = ip + 1;
                    }
                    nTri03++;
                }
            }
        }
        else if (element_type == Z_QUA04)
        {
            forAll(bMesh[patchi], n)
            {
                label nPoints = bMesh[patchi][n].size();
                if (nPoints == 4)
                {
                    for (label i=0; i<nPoints; i++)
                    {
                        label ip = bMesh[patchi][n][i];
                        conn_array[nQuad04][i] = ip + 1;
                    }
                    nQuad04++;
                }
            }

        }
        else if (element_type == Z_NSIDED)
        {
            label nPoly = 0;
            forAll(bMesh[patchi], n)
            {
                label nPoints = bMesh[patchi][n].size();
                if ((nPoints != 3) && (nPoints != 4))
                {
                    conn_array[nPoly++][0] = nPoints;
                }
            }
        }
    }
    else if (part_number == nPatches+2)
    {
        forAll(*sprayPtr, n)
        {
            conn_array[n][0] = n + 1;
        }
    }
    else
    {
        return Z_ERR;
    }

    #ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_part_elements_by_type" << endl;
    #endif

    return Z_OK;
}
