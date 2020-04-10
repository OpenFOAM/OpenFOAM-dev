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
IOdictionary momentumTransport
(
    IOobject
    (
        "momentumTransport",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
const dictionary& MaxwellCoeffs =
    momentumTransport.subDict("laminar").subDict("MaxwellCoeffs");
const scalar nu1 = MaxwellCoeffs.lookup<scalar>("nuM");
const scalar lambda = MaxwellCoeffs.lookup<scalar>("lambda");

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
    fvOptions.subDict("momentumSource").subDict("sources");
const scalar K(vector(gradPDict.subDict("U").lookup("explicit")).x());

dictionary probes(IFstream(runTime.system()/"probes")());
const point location = pointField(probes.lookup("probeLocations"))[0];
