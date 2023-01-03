#include "createControl.H"

int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

scalar convergenceTolerance(stressControl.lookupOrDefault<scalar>("D", 0));

scalar accFac
(
    stressControl.lookupOrDefault<scalar>("accelerationFactor", 1)
);
