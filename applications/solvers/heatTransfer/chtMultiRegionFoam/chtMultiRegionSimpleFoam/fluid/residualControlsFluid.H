// Residual control used
if (residualControlUsed)
{
    bool UConv = false;
    bool p_rghConv = false;
    bool EConv = false;

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

        if (ETol == -1)
        {
            EConv = true;
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

    if (!EConv)
    {
        if (ETol > solvPerfE.initialResidual())
        {
            EConv = true;
        }
    }

    // Check if each field is converged
    if (UConv && p_rghConv && EConv)
    {
        resReachedFluid = true;
    }
}
