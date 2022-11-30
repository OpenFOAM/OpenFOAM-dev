Foam::word regionName;

if (args.optionReadIfPresent("region", regionName))
{
    Foam::Info
        << "Create polyMesh " << regionName << " for time = "
        << runTime.name() << Foam::nl << Foam::endl;
}
else
{
    regionName = Foam::polyMesh::defaultRegion;
    Foam::Info
        << "Create polyMesh for time = "
        << runTime.name() << Foam::nl << Foam::endl;
}

Foam::polyMesh mesh
(
    Foam::IOobject
    (
        regionName,
        runTime.name(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
