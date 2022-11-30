if (nVar >= Num_variables - nSprayVariables)
{
    return Z_UNDEF;
}


typeIOobject<volVectorField> fieldObjectPtr
(
    fieldNames[var2field[nVar]],
    runTime.name(),
    mesh,
    IOobject::NO_READ
);

if (!fieldObjectPtr.headerOk())
{
    return Z_UNDEF;
}

IOobject fieldObject
(
    fieldNames[var2field[nVar]],
    runTime.name(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

volVectorField vectorField
(
    fieldObject,
    mesh
);

const cellShapeList& cellShapes = meshPtr->cellShapes();

// hexa's
if (which_type == Z_HEX08)
{
    const cellModel& hex = *(cellModeller::lookup("hex"));
    // const cellModel& wedge = *(cellModeller::lookup("wedge"));

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const cellShape& cellShape = cellShapes[n];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == hex) // || (cellModel == wedge))
        {
            var_array[counter++] = vectorField[n][component];
        }
    }
}

// penta's
if (which_type == Z_PEN06)
{
    const cellModel& prism = *(cellModeller::lookup("prism"));

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const cellShape& cellShape = cellShapes[n];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == prism)
        {
            var_array[counter++] = vectorField[n][component];
        }
    }
}

// pyramids's
if (which_type == Z_PYR05)
{
    const cellModel& pyr = *(cellModeller::lookup("pyr"));

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const cellShape& cellShape = cellShapes[n];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == pyr)
        {
            var_array[counter++] = vectorField[n][component];
        }
    }
}


// tet's
if (which_type == Z_TET04)
{
    const cellModel& tet = *(cellModeller::lookup("tet"));

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const cellShape& cellShape = cellShapes[n];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == tet)
        {
            var_array[counter++] = vectorField[n][component];
        }
    }
}

if (which_type == Z_NFACED)
{
    const cellList& cells = meshPtr->cells();

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const labelList& points = cellShapes[n];
        label nFacesInCell = cells[n].size();

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
                var_array[counter++] = vectorField[n][component];
            }
        }
        else
        {
            var_array[counter++] = vectorField[n][component];
        }
    }
}
