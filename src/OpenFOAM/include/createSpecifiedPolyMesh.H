#include "setMeshPath.H"
#include "setRegionName.H"

Foam::polyMesh mesh
(
    Foam::IOobject
    (
        regionName,
        runTime.name(),
        meshPath,
        runTime,
        Foam::IOobject::MUST_READ
    )
);
