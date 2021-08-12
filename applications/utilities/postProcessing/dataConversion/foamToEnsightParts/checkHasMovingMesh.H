// check for "points" in all of the result directories

bool hasMovingMesh = false;
if (timeDirs.size() > 1)
{
    hasMovingMesh = true;
    for (label i=0; i < timeDirs.size() && hasMovingMesh; ++i)
    {
        hasMovingMesh = typeIOobject<pointIOField>
        (
            "points",
            timeDirs[i].name(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ
        ).headerOk();
    }
}
