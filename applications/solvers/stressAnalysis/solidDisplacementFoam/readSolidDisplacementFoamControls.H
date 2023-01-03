{
    const dictionary& stressControl =
        mesh.solution().dict().subDict("stressAnalysis");

    nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);
    convergenceTolerance = stressControl.lookupOrDefault<scalar>("D", 0);
    stressControl.lookup("compactNormalStress") >> compactNormalStress;

    accFac = stressControl.lookupOrDefault<scalar>("accelerationFactor", 1);
}
