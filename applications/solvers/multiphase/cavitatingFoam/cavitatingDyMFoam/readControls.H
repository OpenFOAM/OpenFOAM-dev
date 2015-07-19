#include "readTimeControls.H"

maxAcousticCo = readScalar(runTime.controlDict().lookup("maxAcousticCo"));
correctPhi = pimple.dict().lookupOrDefault<Switch>("correctPhi", true);
