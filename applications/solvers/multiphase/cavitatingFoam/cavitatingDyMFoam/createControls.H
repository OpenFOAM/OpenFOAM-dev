#include "createDyMControls.H"

scalar maxAcousticCo
(
    readScalar(runTime.controlDict().lookup("maxAcousticCo"))
);
