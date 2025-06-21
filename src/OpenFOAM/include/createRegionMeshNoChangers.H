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
    false,
    false
);

mesh.postConstruct(false, fvMesh::stitchType::geometric);
