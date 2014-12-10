    const dictionary& potentialFlow =
        mesh.solutionDict().subDict("potentialFlow");

    const int nNonOrthCorr =
        potentialFlow.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
