const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");

const int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

scalar convergenceTolerance(readScalar(stressControl.lookup("D")));
Switch compactNormalStress(stressControl.lookup("compactNormalStress"));
