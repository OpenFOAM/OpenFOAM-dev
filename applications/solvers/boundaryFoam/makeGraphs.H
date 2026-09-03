const volSymmTensorField R
(
    IOobject
    (
        "R",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->R()
);

runTime.write();

setWriter::New(runTime.controlDict().lookup("graphFormat"))->write
(
    runTime.globalPath()
   /functionObjects::writeFile::outputPrefix
   /args.executable()
   /runTime.name(),

    args.executable(),

    coordSet(true, "y", y),

    "Uf", eval(flowDirection & U.primitiveField())(),

    "nu", turbulence->nu()().primitiveField(),
    "nuEff", turbulence->nuEff()().primitiveField(),
    "k", turbulence->k()().primitiveField(),
    "epsilon", turbulence->epsilon()().primitiveField(),

    "Rff", eval(flowDirection & R.primitiveField() & flowDirection)(),
    "Rww", eval(wallNormal & R.primitiveField() & wallNormal)(),
    "Rfw", eval(flowDirection & R.primitiveField() & wallNormal)(),

    "u", eval(sqrt(mag(R.primitiveField().component(symmTensor::XX))))(),
    "v", eval(sqrt(mag(R.primitiveField().component(symmTensor::YY))))(),
    "w", eval(sqrt(mag(R.primitiveField().component(symmTensor::ZZ))))(),

    "uv", eval(R.primitiveField().component(symmTensor::XY))(),

    "gammaDot", eval(mag(fvi::grad(U)))().primitiveField()
);
