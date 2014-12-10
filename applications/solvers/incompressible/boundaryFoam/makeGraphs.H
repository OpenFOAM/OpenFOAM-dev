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
    turbulence->R()
);

runTime.write();

const word& gFormat = runTime.graphFormat();

makeGraph(y, flowDirection & U, "Uf", gFormat);

makeGraph(y, turbulence->nu(), gFormat);
makeGraph(y, turbulence->k(), gFormat);
makeGraph(y, turbulence->epsilon(), gFormat);

makeGraph(y, flowDirection & R & flowDirection, "Rff", gFormat);
makeGraph(y, wallNormal & R & wallNormal, "Rww", gFormat);
makeGraph(y, flowDirection & R & wallNormal, "Rfw", gFormat);

makeGraph(y, sqrt(mag(R.component(symmTensor::XX))), "u", gFormat);
makeGraph(y, sqrt(mag(R.component(symmTensor::YY))), "v", gFormat);
makeGraph(y, sqrt(mag(R.component(symmTensor::ZZ))), "w", gFormat);
makeGraph(y, R.component(symmTensor::XY), "uv", gFormat);

makeGraph(y, mag(fvc::grad(U)), "gammaDot", gFormat);
