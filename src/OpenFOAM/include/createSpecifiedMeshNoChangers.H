#include "setMeshPath.H"
#include "setRegionName.H"

Foam::fvMesh mesh
(
    Foam::IOobject
    (
        regionName,
        runTime.name(),
        meshPath,
        runTime,
        Foam::IOobject::MUST_READ
    ),
    false
);

mesh.postConstruct(false, fvMesh::stitchType::geometric);
