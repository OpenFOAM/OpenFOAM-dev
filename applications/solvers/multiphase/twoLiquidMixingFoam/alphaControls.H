const dictionary& alphaControls = mesh.solution().solverDict(alpha1.name());

label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));
