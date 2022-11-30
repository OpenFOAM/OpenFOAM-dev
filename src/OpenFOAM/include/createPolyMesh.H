Foam::Info
    << "Create polyMesh for time = "
    << runTime.name() << Foam::nl << Foam::endl;

Foam::polyMesh mesh
(
    Foam::IOobject
    (
        Foam::polyMesh::defaultRegion,
        runTime.name(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
