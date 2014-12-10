#include "readTimeControls.H"

scalar maxAcousticCo
(
    readScalar(runTime.controlDict().lookup("maxAcousticCo"))
);

bool correctPhi =
    pimple.dict().lookupOrDefault<Switch>("correctPhi", true);
