#include "readTimeControls.H"

bool correctPhi = pimple.dict().lookupOrDefault<Switch>
(
    "correctPhi",
    !isA<staticFvMesh>(mesh)
);

bool checkMeshCourantNo =
    pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false);
