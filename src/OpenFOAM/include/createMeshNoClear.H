//
// createMeshNoClear.H
// ~~~~~~~~~~~~~~~~~~~
// currently identical to createMesh.H

    Foam::Info
        << "Create mesh, no clear-out for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
