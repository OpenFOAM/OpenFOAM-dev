#include "setRegionName.H"

Foam::fvMesh mesh
(
    Foam::IOobject
    (
        regionName,
        runTime.name(),
        runTime,
        Foam::IOobject::MUST_READ
    ),
    false
);

mesh.postConstruct(false, false, fvMesh::stitchType::geometric);
