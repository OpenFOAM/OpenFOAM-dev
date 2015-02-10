int USERD_get_var_by_component
(
    int which_variable,
    int which_part,
    int var_type,
    int which_type,
    int imag_data,
    int component,
    float *var_array
)
{
#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_var_by_component" << endl
        << "which_variable = " << which_variable << endl
        << "which_part = " << which_part << endl
        << "var_type = " << var_type << endl
        << "which_type = " << which_type << endl
        << "component = " << component << endl
        << flush;
#endif

    label nVar = which_variable - 1;

    Time& runTime = *runTimePtr;

    fvMesh& mesh = *meshPtr;
    const cellShapeList& cells = mesh.cellShapes();

    label nCells = cells.size();

    if (var_type == Z_SCALAR)
    {
        if (which_part == 1)
        {
            #include "getFieldScalar.H"
        }
        else if (which_part < nPatches+2)
        {
            #include "getPatchFieldScalar.H"
        }
        else if (which_part == nPatches+2)
        {
            #include "getLagrangianScalar.H"
        }
        else
        {
            return Z_ERR;
        }
    }
    else if (var_type == Z_VECTOR)
    {
        if (which_part == 1)
        {
            #include "getFieldVector.H"
        }
        else if (which_part < nPatches+2)
        {
            #include "getPatchFieldVector.H"
        }
        else if (which_part == nPatches+2)
        {
            #include "getLagrangianVector.H"
        }
        else
        {
            return Z_ERR;
        }

    }
    else if (var_type == Z_TENSOR9)
    {
        // all tensor are treated as asymmetric tensors here

        if (which_part == 1)
        {
            #include "getFieldTensor.H"
        }
        else if (which_part < nPatches+2)
        {
            #include "getPatchFieldTensor.H"
        }
        else if (which_part == nPatches+2)
        {
            return Z_UNDEF;
        }
        else
        {
            return Z_ERR;
        }

    }
    else
    {
        return Z_UNDEF;
    }

#ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_var_by_component" << endl
        << flush;
#endif

    return Z_OK;
}
