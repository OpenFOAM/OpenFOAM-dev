Foam::word regionName;

if (args.optionReadIfPresent("region", regionName))
{
    Foam::Info
        << "Create mesh " << regionName << " for time = "
        << runTime.name() << Foam::nl << Foam::endl;
}
else
{
    regionName = Foam::fvMesh::defaultRegion;
    Foam::Info
        << "Create mesh for time = "
        << runTime.name() << Foam::nl << Foam::endl;
}

Foam::fvMesh mesh
(
    Foam::IOobject
    (
        regionName,
        runTime.name(),
        runTime,
        Foam::IOobject::MUST_READ
    ),
    false
);
