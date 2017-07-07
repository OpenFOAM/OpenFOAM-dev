{
    IOobject io
    (
        "points",
        runTime.timeName(),
        polyMesh::meshSubDir,
        mesh
    );

    if (io.typeHeaderOk<pointIOField>(true))
    {
        // Read new points
        io.readOpt() = IOobject::MUST_READ;
        mesh.movePoints(pointIOField(io));
    }
}
