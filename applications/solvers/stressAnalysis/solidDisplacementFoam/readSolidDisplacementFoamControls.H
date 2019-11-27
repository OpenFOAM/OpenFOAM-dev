nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);
convergenceTolerance = stressControl.lookup<scalar>("D");
stressControl.lookup("compactNormalStress") >> compactNormalStress;
