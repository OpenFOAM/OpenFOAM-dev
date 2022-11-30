// Currently identical to createMesh.H

Foam::Info
    << "Create mesh, no clear-out for time = "
    << runTime.name() << Foam::nl << Foam::endl;

Foam::fvMesh mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.name(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
