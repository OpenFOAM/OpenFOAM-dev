const dictionary& residualControl =
    mesh.solutionDict().subDict("PIMPLE").subOrEmptyDict("residualControl");

scalar UTol = -1;
scalar heTol = -1;
scalar p_rghTol = -1;

if (!residualControl.empty())
{
    if (!residualControl.subOrEmptyDict("U").empty())
    {
        UTol = readScalar(residualControl.subDict("U").lookup("tolerance"));
    }

    if (!residualControl.subOrEmptyDict("p_rgh").empty())
    {
        p_rghTol =
            readScalar
            (
                residualControl.subDict("p_rgh").lookup("tolerance")
            );
    }

    if (!residualControl.subOrEmptyDict("h").empty())
    {
        heTol = readScalar(residualControl.subDict("h").lookup("tolerance"));
    }

    // Residual control used?
    if (UTol != -1 || heTol != -1 || p_rghTol != -1)
    {
        residualControlUsed = true;
        resReachedFluid = false;
    }

    if (firstIteration)
    {
        firstIteration = false;

        if (UTol != -1 || heTol != -1 || p_rghTol != -1)
        {
            Info<< "PIMPLE convergence criteria\n";

            if (UTol != -1)
            {
                Info<< "   field U      tolerance " << UTol << "\n";
            }
            if (heTol != -1)
            {
                Info<< "   field h      tolerance " << heTol << "\n";
            }
            if (p_rghTol != -1)
            {
                Info<< "   field p_rgh  tolerance " << p_rghTol << "\n";
            }

            Info<< endl;
        }
    }
}
