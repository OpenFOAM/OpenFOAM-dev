const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
