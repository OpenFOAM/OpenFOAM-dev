const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");

scalar accFac(readScalar(stressControl.lookup("accelerationFactor")));
