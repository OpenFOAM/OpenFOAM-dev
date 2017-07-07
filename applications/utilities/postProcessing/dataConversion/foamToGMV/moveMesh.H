IOobject ioPoints
(
    "points",
    runTime.timeName(),
    mesh.name(),
    mesh
);

if (ioPoints.typeHeaderOk<pointIOField>(true))
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

    mesh.movePoints(newPoints);
}
