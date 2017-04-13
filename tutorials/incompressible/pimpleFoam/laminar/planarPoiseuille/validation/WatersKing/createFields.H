Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
const scalar nu2 =
    dimensionedScalar
    (
        "nu",
        dimViscosity,
        transportProperties.lookup("nu")
    ).value();

Info<< "Reading viscoelastic properties\n" << endl;
IOdictionary turbulenceProperties
(
    IOobject
    (
        "turbulenceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
const dictionary& MaxwellCoeffs =
    turbulenceProperties.subDict("laminar").subDict("MaxwellCoeffs");
const scalar nu1 = readScalar(MaxwellCoeffs.lookup("nuM"));
const scalar lambda = readScalar(MaxwellCoeffs.lookup("lambda"));

const scalar rho = 1;

Info<< "Reading pressure gradient\n" << endl;
IOdictionary fvOptions
(
    IOobject
    (
        "fvOptions",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
const dictionary& gradPDict =
    fvOptions.subDict("momentumSource").subDict("injectionRateSuSp");
const scalar K =
    Tuple2<vector, scalar>(gradPDict.lookup("U")).first().x();

dictionary probes(IFstream(runTime.system()/"probes")());
const point location = pointField(probes.lookup("probeLocations"))[0];
