volSymmTensorField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->sigma()
);

runTime.write();

setWriter::New(runTime.graphFormat())->write
(
    runTime.globalPath()
   /functionObjects::writeFile::outputPrefix
   /args.executable()
   /runTime.timeName(),

    args.executable(),

    coordSet(true, "y", y),

    "Uf", (flowDirection & U)().primitiveField(),

    "nu", turbulence->nu()().primitiveField(),
    "k", turbulence->nu()().primitiveField(),
    "epsilon", turbulence->nu()().primitiveField(),

    "Rff", (flowDirection & R & flowDirection)().primitiveField(),
    "Rww", (wallNormal & R & wallNormal)().primitiveField(),
    "Rfw", (flowDirection & R & wallNormal)().primitiveField(),

    "u", sqrt(mag(R.component(symmTensor::XX)))().primitiveField(),
    "v", sqrt(mag(R.component(symmTensor::YY)))().primitiveField(),
    "w", sqrt(mag(R.component(symmTensor::ZZ)))().primitiveField(),

    "uv", R.component(symmTensor::XY)().primitiveField(),

    "gammaDot", mag(fvc::grad(U))().primitiveField()
);
