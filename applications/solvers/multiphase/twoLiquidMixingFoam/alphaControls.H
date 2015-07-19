const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
