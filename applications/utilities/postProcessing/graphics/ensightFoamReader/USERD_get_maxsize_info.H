int USERD_get_maxsize_info
(
    int *max_number_of_nodes,
    int *max_number_of_elements[Z_MAXTYPE],
    int *max_ijk_dimensions[3]
)
{
    return Z_ERR;

#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_maxsize_info" << endl;
#endif

    label maxNPoints = 0;
    label maxNParcels = 0;

    label nPen06Max = 0;
    label nHex08Max = 0;
    label nPyr05Max = 0;
    label nTet04Max = 0;

    Info<< "Checking all time steps for EnSight memory allocation purpose. "
        << "This can take some time." << endl;

    for (label timeI=1; timeI < timeDirs.size(); ++timeI)
    {

        label nPen06 = 0;
        label nHex08 = 0;
        label nPyr05 = 0;
        label nTet04 = 0;

        runTimePtr->setTime(timeDirs[timeI], timeI);

        Info<< "Checking time = " << runTimePtr->timeName() << endl;

        const cellShapeList& cells = meshPtr->cellShapes();

        const label nPoints = meshPtr->nPoints();
        const label nCells  = cells.size();

        maxNPoints = max(maxNPoints, nPoints);

        for (label n=0; n<nCells;n++)
        {
            label nFaces = cells[n].nFaces();
            const labelList& points = cells[n];

            if ((nFaces == 6) && (points.size() == 8))
            {
                nHex08++;
            }
            else if ((nFaces == 5) && (points.size() == 6))
            {
                nPen06++;
            }
            else if ((nFaces == 5) && (points.size() == 5))
            {
                nPyr05++;
            }
            else if ((nFaces == 4) && (points.size() == 4))
            {
                nTet04++;
            }
        }

        nPen06Max = max(nPen06Max, nPen06);
        nHex08Max = max(nHex08Max, nHex08);
        nPyr05Max = max(nPyr05Max, nPyr05);
        nTet04Max = max(nTet04Max, nTet04);

        if (Numparts_available > 1)
        {
            // Get the maximum number of spray parcels
            // and store it
            Cloud<passiveParticle> lagrangian(*meshPtr, cloud::defaultName);

            if (lagrangian.size() > nMaxParcels)
            {
                nMaxParcels = lagrangian.size();
            }
        }
    }

    max_number_of_nodes[0] = maxNPoints;
    max_number_of_elements[0][Z_HEX08] = nHex08Max;
    max_number_of_elements[0][Z_PEN06] = nPen06Max;
    max_number_of_elements[0][Z_PYR05] = nPyr05Max;
    max_number_of_elements[0][Z_TET04] = nTet04Max;

    if (Numparts_available > 1)
    {
        max_number_of_nodes[1] = maxNParcels;
        max_number_of_elements[1][Z_POINT] = maxNParcels;
    }

#ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_maxsize_info" << endl;
#endif

    return Z_OK;
}
