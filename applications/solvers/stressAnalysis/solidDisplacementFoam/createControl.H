const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");

Switch compactNormalStress(stressControl.lookup("compactNormalStress"));
