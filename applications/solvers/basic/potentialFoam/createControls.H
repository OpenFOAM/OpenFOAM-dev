const dictionary& potentialFlow
(
    mesh.solution().dict().subDict("potentialFlow")
);

const int nNonOrthCorr
(
    potentialFlow.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0)
);
