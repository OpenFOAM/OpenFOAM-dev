if (residualControlUsed)
{
    bool UConv = false;
    bool p_rghConv = false;
    bool heConv = false;

    // Check which field is not used for control
    {
        if (UTol == -1 || !momentumPredictor)
        {
            UConv =  true;
        }

        if (p_rghTol == -1)
        {
            p_rghConv = true;
        }

        if (heTol == -1)
        {
            heConv = true;
        }
    }

    // Get the last initial residual of the solvers
    if (momentumPredictor && !UConv)
    {
        if (UTol > cmptMax(solvPerfU.initialResidual()))
        {
            UConv = true;
        }
    }

    if (!p_rghConv)
    {
        if (p_rghTol > solvPerfp_rgh.initialResidual())
        {
            p_rghConv = true;
        }
    }

    if (!heConv)
    {
        if (heTol > solvPerfhe.initialResidual())
        {
            heConv = true;
        }
    }

    // Check if each field is converged
    if (UConv && p_rghConv && heConv)
    {
        resReachedFluid = true;
    }
}
