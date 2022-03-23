const dictionary& stressControl =
    mesh.solution().dict().subDict("stressAnalysis");

Switch compactNormalStress(stressControl.lookup("compactNormalStress"));
