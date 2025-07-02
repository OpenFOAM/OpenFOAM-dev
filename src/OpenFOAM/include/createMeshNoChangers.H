Foam::Info
    << "Create mesh for time = "
    << runTime.name() << Foam::nl << Foam::endl;

Foam::fvMesh mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.name(),
        runTime,
        Foam::IOobject::MUST_READ
    ),
    false
);

mesh.postConstruct(false, false, fvMesh::stitchType::geometric);
