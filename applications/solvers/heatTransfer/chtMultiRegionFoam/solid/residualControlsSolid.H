// Residual control used
if (residualControlUsed)
{
    if (hTol > solvPerfh.initialResidual())
    {
        resReachedSolid = true;
    }
}
