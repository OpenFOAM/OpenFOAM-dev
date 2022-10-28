Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiuMulticomponentThermo> pThermo
(
    psiuMulticomponentThermo::New(mesh)
);
psiuMulticomponentThermo& thermo = pThermo();
thermo.validate(args.executable(), "ha", "ea");

basicCombustionMixture& composition = thermo.composition();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.renameRho()
);

volScalarField& p = thermo.p();

volScalarField& b = composition.Y("b");
Info<< "min(b) = " << min(b).value() << endl;


Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

mesh.schemes().setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermoThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermoThermophysicalTransportModel::New(turbulence(), thermo)
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

Info<< "Creating field Xi\n" << endl;
volScalarField Xi
(
    IOobject
    (
        "Xi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Creating the unstrained laminar flame speed\n" << endl;
autoPtr<laminarFlameSpeed> unstrainedLaminarFlameSpeed
(
    laminarFlameSpeed::New(thermo)
);


Info<< "Reading strained laminar flame speed field Su\n" << endl;
volScalarField Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

dimensionedScalar SuMin = 0.01*Su.average();
dimensionedScalar SuMax = 4*Su.average();

Info<< "Calculating turbulent flame speed field St\n" << endl;
volScalarField St
(
    IOobject
    (
        "St",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Xi*Su
);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

if (composition.contains("ft"))
{
    fields.add(composition.Y("ft"));
}

fields.add(b);
fields.add(thermo.he());
fields.add(thermo.heu());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
