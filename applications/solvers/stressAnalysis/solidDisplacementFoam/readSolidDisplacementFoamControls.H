nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);
convergenceTolerance = readScalar(stressControl.lookup("D"));
stressControl.lookup("compactNormalStress") >> compactNormalStress;
