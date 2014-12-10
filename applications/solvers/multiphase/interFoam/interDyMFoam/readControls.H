#include "readTimeControls.H"

bool correctPhi
(
    pimple.dict().lookupOrDefault<Switch>("correctPhi", true)
);

bool checkMeshCourantNo
(
    pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false)
);

bool moveMeshOuterCorrectors
(
    pimple.dict().lookupOrDefault<Switch>("moveMeshOuterCorrectors", false)
);
