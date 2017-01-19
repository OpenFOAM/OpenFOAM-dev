//======================================================================
int USERD_get_var_value_at_specific
(
    int which_var,
    int which_node_or_elem,
    int which_part,
    int which_elem_type,
    int time_step,
    float values[3],
    int imag_data
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_var_value_at_specific" << endl
        << flush;
#endif
    // Not sure if it is 0 or 1 based

    label nNode = which_node_or_elem;
    label nVar = which_var - 1;

    fvMesh& mesh = *meshPtr;

    if (nVar < Num_variables - nSprayVariables)
    {
        Time& runTime = *runTimePtr;

        IOobject fieldObject
        (
            fieldNames[var2field[nVar]],
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (isScalar[nVar])
        {
            volScalarField scalarField(fieldObject,mesh);
            values[0] = scalarField[nNode];
        }
        else if (isVector[nVar])
        {
            volVectorField vectorField(fieldObject,mesh);
            values[0] = vectorField[nNode].x();
            values[1] = vectorField[nNode].y();
            values[2] = vectorField[nNode].z();
        }
        else
        {
            Info<< "ERROR in USERD_get_variable_value_at_specific. "
                << "No available variable???"
                << endl;
            return Z_ERR;
        }
    }
    else
    {
        Info<< "This functionality is not implemented yet."
            << endl;
        return Z_ERR;
    }


#ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_var_value_at_specific" << endl
        << flush;
#endif
    return Z_OK;
}
