const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));
