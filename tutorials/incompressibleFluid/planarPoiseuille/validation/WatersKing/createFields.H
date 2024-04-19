Info<< "Reading physicalProperties\n" << endl;

autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));

const scalar nu2 =
    dimensionedScalar
    (
        "nu",
        dimKinematicViscosity,
        viscosity->lookup("nu")
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
IOdictionary fvModels
(
    IOobject
    (
        "fvModels",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
const dictionary& gradPDict =
    fvModels.subDict("momentumSource").subDict("sources");
const scalar K(vector(gradPDict.subDict("U").lookup("explicit")).x());

const dictionary& probesDict = runTime.functionObjects().subDict("probes");
const point location = pointField(probesDict.lookup("probeLocations"))[0];
