#include "createControl.H"

int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

scalar convergenceTolerance(readScalar(stressControl.lookup("D")));
