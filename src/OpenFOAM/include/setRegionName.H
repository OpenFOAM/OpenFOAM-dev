Foam::word regionName;

if (args.optionReadIfPresent("region", regionName))
{
    Foam::Info
        << "Create mesh " << regionName << " for time = "
        << runTime.name() << Foam::nl << Foam::endl;
}
else
{
    regionName = Foam::polyMesh::defaultRegion;
    Foam::Info
        << "Create mesh for time = "
        << runTime.name() << Foam::nl << Foam::endl;
}
