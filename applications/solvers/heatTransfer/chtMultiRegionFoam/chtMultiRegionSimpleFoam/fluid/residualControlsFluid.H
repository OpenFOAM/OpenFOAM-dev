// Residual control used
if (residualControlUsed)
{
    bool UConv = false;
    bool p_rghConv = false;
    bool hConv = false;

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

        if (hTol == -1)
        {
            hConv = true;
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

    if (!hConv)
    {
        if (hTol > solvPerfh.initialResidual())
        {
            hConv = true;
        }
    }

    // Check if each field is converged
    if (UConv && p_rghConv && hConv)
    {
        resReachedFluid = true;
    }
}
