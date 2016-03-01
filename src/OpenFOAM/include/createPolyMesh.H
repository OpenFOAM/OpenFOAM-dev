Foam::Info
    << "Create polyMesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

Foam::polyMesh mesh
(
    Foam::IOobject
    (
        Foam::polyMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
