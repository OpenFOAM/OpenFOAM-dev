#include "readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault
(
    "correctPhi",
    correctPhi
);

checkMeshCourantNo = pimple.dict().lookupOrDefault
(
    "checkMeshCourantNo",
    checkMeshCourantNo
);

moveMeshOuterCorrectors = pimple.dict().lookupOrDefault
(
    "moveMeshOuterCorrectors",
    moveMeshOuterCorrectors
);
