// Residual control used?
bool resControlUsed = false;
int nFluidControlled = fluidRegions.size();
int nSolidControlled = solidRegions.size();

// Check wheater there is a single regions that uses residual control
forAll(fluidRegions, i)
{
    if (residualControlUsedFluid[i])
    {
        resControlUsed = true;
        break;
    }
}

forAll(solidRegions, i)
{
    if(residualControlUsedSolid[i])
    {
        resControlUsed = true;
        break;
    }
}

if (resControlUsed)
{
    int nFluidConv = 0;
    int nSolidConv = 0;

    // Sum of all converged regions (Note: if no residual control is used
    // the residualReached* flag is already set to true)
    forAll(fluidRegions, i)
    {
        if (residualReachedFluid[i])
        {
            nFluidConv++;
        }
    }

    forAll(solidRegions, i)
    {
        if (residualReachedSolid[i])
        {
            nSolidConv++;
        }
    }

    if (nFluidConv == nFluidControlled && nSolidConv == nSolidControlled)
    {
        // Activate flag to go to the 'Final' loop using the 'Final'
        // relaxation factors
        allRegionsConverged = true;
    }
}

if (!runTime.run())
{
    Info<< "\nRegions not converged after " << runTime.timeName()
        << " iterations" << endl;
}
else if (runTime.run() && resControlUsed && allRegionsConverged)
{
    Info<< "\nRegions converged after " << runTime.timeName()
        << " iterations" << endl;

    // Save converged time step and end the run
    runTime.writeAndEnd();
}
