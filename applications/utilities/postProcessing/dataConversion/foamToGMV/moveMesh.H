typeIOobject<pointIOField> ioPoints
(
    "points",
    runTime.timeName(),
    mesh.name(),
    mesh
);

if (ioPoints.headerOk())
{
    Info<< "new points available" << endl;
    // Reading new points
    pointIOField newPoints
    (
        IOobject
        (
            "points",
            runTime.timeName(),
            mesh.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    mesh.setPoints(newPoints);
}
